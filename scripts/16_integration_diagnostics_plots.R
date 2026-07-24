#!/usr/bin/env Rscript
# =============================================================================
# 16_integration_analysis.R — Intégration RNA-seq x ChIP-seq (FRDA)   [v2]
# -----------------------------------------------------------------------------
# PRÉREQUIS : 04b_recompute_gene_level.R (correctif du logFC ChIP).
#
# CORRECTIFS v2
#   - dplyr::select() qualifié partout (org.Hs.eg.db masque select)
#   - déduplication RNA par ENTREZID -> supprime le many-to-many
#   - rapporte l'ODDS RATIO et l'EXCÈS (obs - null), pas seulement le % brut :
#     le null n'est PAS à 50% quand les marges sont asymétriques.
#   - écarte les logFC RNA aberrants (|logFC| > 15 = comptage nul dans un groupe)
#
# Usage : export PROJECT_ROOT=$(pwd) && Rscript 16_integration_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(tibble)
  library(ggplot2); library(ggrepel); library(stringr)
  library(org.Hs.eg.db); library(AnnotationDbi)
})

save_plot <- function(p, filename, width=11, height=8) {
  ggsave(file.path(PLOT_DIR, filename), p, width = width, height = height)
}

set.seed(42)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
ANN_DIR  <- file.path(PROJECT_ROOT, "results", "annotation")
INT_DIR  <- file.path(PROJECT_ROOT, "results", "integration")
PLOT_DIR <- file.path(INT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

RNA_FILE <- "/home/apuhart/Downloads/Bulk-RNAseq/results_test/fibro/control_vs_FA/fibro_control_vs_FA_DESeq2_results.csv"

MARK_TYPE <- c(H3K9me3 = "repressive", H3K27me3 = "repressive",
               H3K27ac = "active",     PolII    = "active")

N_PERM      <- 1000
FDR_RNA_TH  <- 0.05
LFC_RNA_TH  <- 1
LFC_RNA_MAX <- 15

PAL <- c(repressive = "#4E79A7", active = "#E15759")

# ---------------------------------------------------------------------------
# 0a. RNA-seq
# ---------------------------------------------------------------------------
cat("\n=========================== CHARGEMENT ===========================\n")

rna_raw <- read_csv(RNA_FILE, show_col_types = FALSE)
names(rna_raw)[1] <- "GENE_ID"

need <- c("log2FoldChange", "pvalue", "padj")
miss <- setdiff(need, names(rna_raw))
if (length(miss)) stop("Colonnes absentes du CSV RNA : ", paste(miss, collapse = ", "))

rna <- rna_raw %>%
  transmute(GENE_ID    = sub("\\..*$", "", as.character(GENE_ID)),
            baseMean   = if ("baseMean" %in% names(rna_raw)) as.numeric(baseMean) else NA_real_,
            logFC_RNA  = as.numeric(log2FoldChange),
            FDR_RNA    = as.numeric(padj),
            PValue_RNA = as.numeric(pvalue)) %>%
  filter(!is.na(GENE_ID), GENE_ID != "")

# |logFC| > 15 -> comptage nul dans un groupe : 2^24 = 16 millions de fois.
# Ce n'est pas un effet biologique, c'est une division par (presque) zéro.
n_extreme <- sum(abs(rna$logFC_RNA) > LFC_RNA_MAX, na.rm = TRUE)
if (n_extreme) {
  cat(sprintf("!! %d gènes avec |logFC_RNA| > %d -> écartés (comptage nul).\n",
              n_extreme, LFC_RNA_MAX))
  rna <- rna %>% filter(is.na(logFC_RNA) | abs(logFC_RNA) <= LFC_RNA_MAX)
}

is_ens <- grepl("^ENSG[0-9]+$", rna$GENE_ID)

safe_select <- function(keys_in, keytype, columns) {
  ok <- intersect(unique(keys_in), AnnotationDbi::keys(org.Hs.eg.db, keytype))
  if (!length(ok)) return(NULL)
  suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = ok,
                                         keytype = keytype, columns = columns))
}
m_sym <- safe_select(rna$GENE_ID[!is_ens], "SYMBOL", c("ENTREZID", "ENSEMBL"))
m_sym <- if (is.null(m_sym)) tibble() else m_sym %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% transmute(GENE_ID = SYMBOL, ENTREZID, SYMBOL)
m_ens <- safe_select(rna$GENE_ID[is_ens], "ENSEMBL", c("ENTREZID", "SYMBOL"))
m_ens <- if (is.null(m_ens)) tibble() else m_ens %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>% transmute(GENE_ID = ENSEMBL, ENTREZID, SYMBOL)

rna <- rna %>%
  left_join(bind_rows(m_sym, m_ens) %>% distinct(GENE_ID, .keep_all = TRUE),
            by = "GENE_ID") %>%
  filter(!is.na(ENTREZID)) %>%
  mutate(ENTREZID = as.character(ENTREZID))

# DÉDUPLICATION : plusieurs ENSG -> même Entrez. Sans ça, l'inner_join avec le
# ChIP crée un many-to-many et duplique des gènes -> comptages faussés partout.
n_before <- nrow(rna)
rna <- rna %>% arrange(FDR_RNA, PValue_RNA) %>% distinct(ENTREZID, .keep_all = TRUE)
cat(sprintf("RNA  : %d gènes (%d doublons Entrez écartés) | %d à FDR<%.2f\n",
            nrow(rna), n_before - nrow(rna),
            sum(rna$FDR_RNA < FDR_RNA_TH, na.rm = TRUE), FDR_RNA_TH))

# ---------------------------------------------------------------------------
# 0b. ChIP + jointure
# ---------------------------------------------------------------------------
files <- list.files(ANN_DIR, pattern = "_gene_level_fixed\\.csv$", full.names = TRUE)
if (!length(files))
  stop("Aucun *_gene_level_fixed.csv dans ", ANN_DIR,
       "\n-> Lance d'abord 04b_recompute_gene_level.R")

chip <- lapply(files, function(f) {
  tag <- sub("_gene_level_fixed$", "", tools::file_path_sans_ext(basename(f)))
  read_csv(f, show_col_types = FALSE) %>%
    mutate(mark = sub("__.*$", "", tag), ENTREZID = as.character(GeneID))
}) %>% bind_rows() %>%
  distinct(mark, ENTREZID, region_class, .keep_all = TRUE)

dat_all <- chip %>%
  inner_join(rna, by = "ENTREZID", relationship = "many-to-one") %>%
  transmute(mark, type = unname(MARK_TYPE[mark]), region_class,
            gene = coalesce(GeneSymbol, SYMBOL),
            ENTREZID, n_regions, mean_signal,
            logFC_chip = logFC, FDR_chip = FDR,
            baseMean, logFC_RNA, PValue_RNA, FDR_RNA) %>%
  filter(!is.na(logFC_chip), !is.na(logFC_RNA))

dat <- dat_all %>% filter(region_class == "All")

# ---------------------------------------------------------------------------
# 0c. DIAGNOSTIC
# ---------------------------------------------------------------------------
cat("\n=========================== DIAGNOSTIC ===========================\n")
diag <- dat %>%
  group_by(mark, type) %>%
  summarise(n_genes = n(),
            max_abs_lfc_chip = round(max(abs(logFC_chip), na.rm = TRUE), 3),
            sig_chip_05 = sum(FDR_chip < 0.05, na.rm = TRUE),
            max_abs_lfc_rna = round(max(abs(logFC_RNA), na.rm = TRUE), 2),
            sig_rna_05 = sum(FDR_RNA < 0.05, na.rm = TRUE),
            .groups = "drop")
print(as.data.frame(diag), row.names = FALSE)

if (any(diag$max_abs_lfc_chip < 0.1))
  stop("\n!! logFC ChIP TOUJOURS ÉCRASÉ : relance 04b_recompute_gene_level.R\n")
cat("\nOK : logFC ChIP dans une plage réaliste.\n")

pct_concordant <- function(lfc_chip, lfc_rna, type) {
  keep <- lfc_chip != 0 & lfc_rna != 0 & !is.na(lfc_chip) & !is.na(lfc_rna)
  s <- sign(lfc_chip[keep]) * sign(lfc_rna[keep])
  if (!length(s)) return(NA_real_)
  100 * mean(if (type == "repressive") s < 0 else s > 0)
}

# ===========================================================================
# 1. CONCORDANCE GLOBALE
# ===========================================================================
# ATTENTION : le % brut dépend des marges. Si 70% des logFC_chip sont négatifs,
# le % concordant attendu par hasard n'est PAS 50%. L'ODDS RATIO est, lui,
# indépendant des marges. Quand les deux divergent, croire l'OR.
# Sens attendu : OR > 1 pour une marque ACTIVE, OR < 1 pour une RÉPRESSIVE.
cat("\n=================== 1. CONCORDANCE GLOBALE ===================\n")

conc <- lapply(split(dat, dat$mark), function(d) {
  d <- d %>% filter(logFC_chip != 0, logFC_RNA != 0)
  tab <- table(chip = sign(d$logFC_chip), rna = sign(d$logFC_RNA))
  if (!all(dim(tab) == c(2, 2)))
    return(tibble(mark = d$mark[1], type = d$type[1], n = nrow(d),
                  OR = NA_real_, fisher_p = NA_real_, pct = NA_real_))
  ft <- fisher.test(tab)
  tibble(mark = d$mark[1], type = d$type[1], n = sum(tab),
         OR = round(unname(ft$estimate), 3),
         fisher_p = ft$p.value,
         pct = round(pct_concordant(d$logFC_chip, d$logFC_RNA, d$type[1]), 1))
}) %>% bind_rows() %>% arrange(desc(pct))

print(as.data.frame(conc %>% mutate(fisher_p = signif(fisher_p, 3))), row.names = FALSE)
cat("\nOR attendu : >1 si marque ACTIVE, <1 si RÉPRESSIVE.\n")

# ===========================================================================
# 2. NULL PAR PERMUTATION
# ===========================================================================
cat(sprintf("\n=================== 2. NULL PAR PERMUTATION (%d) ===================\n", N_PERM))

perm <- lapply(split(dat, dat$mark), function(d) {
  d <- d %>% filter(logFC_chip != 0, logFC_RNA != 0)
  ty  <- d$type[1]
  obs <- pct_concordant(d$logFC_chip, d$logFC_RNA, ty)
  nv  <- replicate(N_PERM, pct_concordant(sample(d$logFC_chip), d$logFC_RNA, ty))
  tibble(mark = d$mark[1], type = ty,
         observed  = round(obs, 2),
         null_mean = round(mean(nv), 2),
         excess    = round(obs - mean(nv), 2),
         null_sd   = round(sd(nv), 3),
         z         = round((obs - mean(nv)) / sd(nv), 1),
         p_perm    = (1 + sum(abs(nv - mean(nv)) >= abs(obs - mean(nv)))) / (N_PERM + 1))
}) %>% bind_rows() %>% arrange(desc(excess))

print(as.data.frame(perm), row.names = FALSE)
write_csv(perm, file.path(INT_DIR, "permutation_null.csv"))

cat("\nLECTURE — deux chiffres, pas un :\n")
cat("  excess = observed - null_mean  -> TAILLE de l'effet (points de %)\n")
cat("  z                              -> SIGNIFICATIVITÉ\n")
cat("Un z élevé + un excess minuscule = effet réel mais négligeable.\n")
cat("Le null n'est PAS à 50% : le comparer à 50% donne des conclusions fausses.\n")

perm_plot <- lapply(split(dat, dat$mark), function(d) {
  d <- d %>% filter(logFC_chip != 0, logFC_RNA != 0)
  ty <- d$type[1]
  tibble(mark = d$mark[1], type = ty,
         null = replicate(300, pct_concordant(sample(d$logFC_chip), d$logFC_RNA, ty)))
}) %>% bind_rows() %>%
  left_join(perm %>% dplyr::select(mark, observed), by = "mark")

p_perm <- ggplot(perm_plot, aes(null)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "white") +
  geom_vline(aes(xintercept = observed, colour = type), linewidth = 1) +
  geom_text(data = perm, aes(x = -Inf, y = Inf,
                             label = sprintf("obs = %.1f%%\nnull = %.1f%%\nexcès = %+.1f\nz = %.1f",
                                             observed, null_mean, excess, z)),
            hjust = -0.1, vjust = 1.2, size = 2.8, inherit.aes = FALSE) +
  scale_colour_manual(values = PAL) +
  facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free") +
  labs(x = "% concordant sous permutation", y = "Fréquence",
       title = "Concordance observée vs null par permutation",
       subtitle = "Barre = observé | histogramme = null (étiquettes mélangées)") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(PLOT_DIR, "02_permutation_null.pdf"), p_perm, width = 11, height = 8)

# ===========================================================================
# 3. STRATIFICATION PAR ABONDANCE
# ===========================================================================
cat("\n=================== 3. STRATIFICATION PAR ABONDANCE ===================\n")

strat <- dat %>%
  group_by(mark) %>% mutate(q_signal = ntile(mean_signal, 5)) %>%
  group_by(mark, type, q_signal) %>% filter(n() >= 30) %>%
  summarise(n = n(),
            pct = round(pct_concordant(logFC_chip, logFC_RNA, type[1]), 1),
            .groups = "drop")
print(as.data.frame(strat), row.names = FALSE)
write_csv(strat, file.path(INT_DIR, "concordance_by_abundance.csv"))

p_strat <- ggplot(strat, aes(factor(q_signal), pct, fill = type)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey30") +
  geom_text(aes(label = paste0(pct, "%")), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = PAL) +
  facet_wrap(~ paste0(mark, " (", type, ")")) +
  coord_cartesian(ylim = c(0, 80)) +
  labs(x = "Quintile de signal ChIP (1 = faible couverture → 5 = forte)",
       y = "% concordant",
       title = "La concordance tient-elle à toutes les abondances ?",
       subtitle = "Si l'effet ne vit que dans certains quintiles = artefact") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(PLOT_DIR, "03_concordance_by_abundance.pdf"), p_strat, width = 10, height = 7)

# ===========================================================================
# 4. PROMOTEUR vs DISTAL
# ===========================================================================
cat("\n=================== 4. PROMOTEUR vs DISTAL ===================\n")

by_reg <- dat_all %>%
  filter(region_class %in% c("Promoter", "Genic", "Distal")) %>%
  group_by(mark, type, region_class) %>% filter(n() >= 50) %>%
  summarise(n = n(),
            pct = round(pct_concordant(logFC_chip, logFC_RNA, type[1]), 1),
            mean_lfc_chip = round(mean(logFC_chip, na.rm = TRUE), 3),
            .groups = "drop")

if (nrow(by_reg)) {
  print(as.data.frame(by_reg), row.names = FALSE)
  write_csv(by_reg, file.path(INT_DIR, "concordance_by_region_class.csv"))
  
  p_reg <- ggplot(by_reg, aes(region_class, pct, fill = type)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 50, linetype = 2, colour = "grey30") +
    geom_text(aes(label = paste0(pct, "%\nn=", n)), vjust = -0.3, size = 2.6) +
    scale_fill_manual(values = PAL) +
    facet_wrap(~ paste0(mark, " (", type, ")")) +
    coord_cartesian(ylim = c(0, 85)) +
    labs(x = NULL, y = "% concordant",
         title = "Concordance par classe de région",
         subtitle = "Promoteur = régulation directe | Distal = paysage enhancer") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(PLOT_DIR, "04_concordance_by_region.pdf"), p_reg, width = 10, height = 7)
} else cat("Pas assez de régions classées.\n")

# ===========================================================================
# 5. TENDANCE BINNÉE
# ===========================================================================
binned <- dat %>%
  group_by(mark) %>% mutate(bin = ntile(logFC_RNA, 10)) %>%
  group_by(mark, type, bin) %>%
  summarise(mean_chip = mean(logFC_chip, na.rm = TRUE),
            se_chip = sd(logFC_chip, na.rm = TRUE) / sqrt(n()),
            n = n(), .groups = "drop")

p_bin <- ggplot(binned, aes(bin, mean_chip)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_errorbar(aes(ymin = mean_chip - 1.96 * se_chip,
                    ymax = mean_chip + 1.96 * se_chip), width = 0.2, colour = "grey50") +
  geom_line(aes(colour = type), linewidth = 0.8) +
  geom_point(aes(colour = type), size = 2.5) +
  scale_x_continuous(breaks = 1:10) +
  scale_colour_manual(values = PAL) +
  facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_y") +
  labs(x = "Décile de log2FC RNA (1 = plus réprimé → 10 = plus induit)",
       y = "log2FC ChIP moyen (± IC95)",
       title = "Le signal chromatinien suit-il l'expression ?",
       subtitle = "Répressive : pente descendante attendue | Active : pente montante") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(PLOT_DIR, "05_binned_trend.pdf"), p_bin, width = 11, height = 8)

# ===========================================================================
# 6. SCATTER
# ===========================================================================
lab <- dat %>% group_by(mark) %>%
  filter(gene == "FXN" | rank(FDR_RNA, ties.method = "first") <= 8) %>% ungroup()

p_sc <- ggplot(dat, aes(logFC_chip, logFC_RNA)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(colour = FDR_RNA < FDR_RNA_TH), alpha = 0.35, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "#E15759", linewidth = 0.7) +
  geom_text_repel(data = lab, aes(label = gene), size = 2.4,
                  max.overlaps = 20, min.segment.length = 0) +
  geom_point(data = filter(dat, gene == "FXN"),
             colour = "black", size = 3, shape = 21, stroke = 1.1, fill = NA) +
  scale_colour_manual(values = c(`FALSE` = "grey78", `TRUE` = "#4E79A7"),
                      name = paste0("RNA FDR < ", FDR_RNA_TH)) +
  facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_x") +
  labs(x = "log2FC ChIP (FRDA vs Control)", y = "log2FC RNA (FRDA vs Control)",
       title = "Chromatine vs expression", subtitle = "Cercle noir = FXN") +
  theme_bw(base_size = 11)
ggsave(file.path(PLOT_DIR, "06_scatter.pdf"), p_sc, width = 11, height = 8)

# ===========================================================================
# 7. BOXPLOT + Kruskal-Wallis
# ===========================================================================
dat_box <- dat %>%
  mutate(RNA_class = factor(case_when(
    FDR_RNA < FDR_RNA_TH & logFC_RNA > 0 ~ "Up",
    FDR_RNA < FDR_RNA_TH & logFC_RNA < 0 ~ "Down",
    TRUE ~ "Stable"), levels = c("Down", "Stable", "Up")))

n_lab <- dat_box %>% count(mark, type, RNA_class)

p_box <- ggplot(dat_box, aes(RNA_class, logFC_chip, fill = RNA_class)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_boxplot(outlier.size = 0.2, notch = TRUE) +
  geom_text(data = n_lab, aes(x = RNA_class, y = -Inf, label = paste0("n=", n)),
            vjust = -0.6, size = 2.5, inherit.aes = FALSE) +
  facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_y") +
  scale_fill_manual(values = c(Down = "#4E79A7", Stable = "grey85", Up = "#E15759")) +
  labs(x = "Catégorie RNA", y = "log2FC ChIP",
       title = "Changement chromatinien selon l'état transcriptionnel",
       subtitle = "Encoches disjointes = médianes significativement différentes") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(PLOT_DIR, "07_boxplot_by_RNA_class.pdf"), p_box, width = 10, height = 7)

cat("\n=================== 7. TEST Down/Stable/Up ===================\n")
kw <- lapply(split(dat_box, dat_box$mark), function(d) {
  if (n_distinct(d$RNA_class) < 2) return(NULL)
  k <- kruskal.test(logFC_chip ~ RNA_class, data = d)
  tibble(mark = d$mark[1], type = d$type[1],
         chi2 = round(unname(k$statistic), 1), p = signif(k$p.value, 3),
         med_down   = round(median(d$logFC_chip[d$RNA_class == "Down"], na.rm = TRUE), 3),
         med_stable = round(median(d$logFC_chip[d$RNA_class == "Stable"], na.rm = TRUE), 3),
         med_up     = round(median(d$logFC_chip[d$RNA_class == "Up"], na.rm = TRUE), 3),
         n_down = sum(d$RNA_class == "Down"), n_up = sum(d$RNA_class == "Up"))
}) %>% bind_rows()
if (nrow(kw)) print(as.data.frame(kw), row.names = FALSE)

# ===========================================================================
# 8. ROC SIGNÉE
# ===========================================================================
auc_tab <- tibble()
if (requireNamespace("pROC", quietly = TRUE)) {
  suppressPackageStartupMessages(library(pROC))
  cat("\n=================== 8. ROC SIGNÉE (direction) ===================\n")
  roc_dir <- file.path(PLOT_DIR, "ROC"); dir.create(roc_dir, showWarnings = FALSE)
  
  auc_tab <- lapply(split(dat, dat$mark), function(d) {
    dd <- d %>% filter(FDR_RNA < FDR_RNA_TH, abs(logFC_RNA) > LFC_RNA_TH)
    if (nrow(dd) < 20 || n_distinct(sign(dd$logFC_RNA)) < 2) {
      cat(sprintf("  [%s] n=%d gènes RNA-sig -> trop peu pour une ROC.\n",
                  d$mark[1], nrow(dd)))
      return(NULL)
    }
    resp <- factor(ifelse(dd$logFC_RNA > 0, "Up", "Down"), levels = c("Down", "Up"))
    pred <- if (dd$type[1] == "repressive") -dd$logFC_chip else dd$logFC_chip
    r <- pROC::roc(resp, pred, quiet = TRUE, direction = "<")
    pdf(file.path(roc_dir, paste0("ROC_", dd$mark[1], ".pdf")), width = 5, height = 5)
    plot(r, legacy.axes = TRUE, col = PAL[[dd$type[1]]], lwd = 3,
         main = sprintf("%s (%s)\nAUC = %.3f  |  n = %d",
                        dd$mark[1], dd$type[1], as.numeric(pROC::auc(r)), nrow(dd)))
    abline(0, 1, lty = 2, col = "grey60"); dev.off()
    ci <- as.numeric(pROC::ci.auc(r))
    tibble(mark = dd$mark[1], type = dd$type[1], n = nrow(dd),
           AUC = round(as.numeric(pROC::auc(r)), 3),
           CI_low = round(ci[1], 3), CI_high = round(ci[3], 3))
  }) %>% bind_rows()
  
  if (nrow(auc_tab)) {
    print(as.data.frame(auc_tab), row.names = FALSE)
    write_csv(auc_tab, file.path(INT_DIR, "AUC_signed.csv"))
    cat("\nSi l'IC95 recouvre 0.5, la marque ne prédit pas le sens.\n")
    
    p_auc <- ggplot(auc_tab, aes(reorder(mark, AUC), AUC, fill = type)) +
      geom_col(width = 0.6) +
      geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15) +
      geom_hline(yintercept = 0.5, linetype = 2) +
      geom_text(aes(label = paste0(AUC, "\nn=", n)), hjust = -0.2, size = 3) +
      scale_fill_manual(values = PAL) +
      coord_flip(ylim = c(0.3, 1)) +
      labs(x = NULL, y = "ROC AUC (± IC95)",
           title = "La chromatine prédit-elle le SENS du changement RNA ?") +
      theme_bw(base_size = 11)
    ggsave(file.path(PLOT_DIR, "08_AUC_signed.pdf"), p_auc, width = 7, height = 4.5)
  }
} else cat("\n[skip] package pROC absent.\n")

# ===========================================================================
# 9. PARTITION DIRECT / INDIRECT + GO
# ===========================================================================
cat("\n=================== 9. DIRECT / INDIRECT ===================\n")

FOCUS_MARK <- "H3K27ac"
d27 <- dat %>% filter(mark == FOCUS_MARK)

if (nrow(d27)) {
  chip_thr <- quantile(abs(d27$logFC_chip), 0.75, na.rm = TRUE)
  cat(sprintf("Seuil |logFC_chip| (Q75) = %.3f\n", chip_thr))
  
  d27 <- d27 %>%
    mutate(rna_sig  = !is.na(FDR_RNA) & FDR_RNA < FDR_RNA_TH & abs(logFC_RNA) > LFC_RNA_TH,
           chip_chg = abs(logFC_chip) >= chip_thr,
           group = case_when(
             rna_sig & logFC_RNA < 0 & chip_chg & logFC_chip < 0 ~ "a_direct_down",
             rna_sig & logFC_RNA > 0 & chip_chg & logFC_chip > 0 ~ "a_direct_up",
             rna_sig & !chip_chg                                  ~ "b_indirect",
             !rna_sig & chip_chg                                  ~ "c_chip_only",
             rna_sig & chip_chg                                   ~ "d_discordant",
             TRUE ~ "e_ns"))
  
  print(as.data.frame(count(d27, group)), row.names = FALSE)
  write_csv(d27 %>% dplyr::select(gene, ENTREZID, logFC_RNA, FDR_RNA,
                                  logFC_chip, FDR_chip, group),
            file.path(INT_DIR, sprintf("%s_direct_indirect_groups.csv", FOCUS_MARK)))
  
  if (requireNamespace("clusterProfiler", quietly = TRUE)) {
    suppressPackageStartupMessages(library(clusterProfiler))
    universe_e <- unique(d27$ENTREZID)
    go_res <- list()
    for (g in c("a_direct_down", "a_direct_up", "b_indirect", "c_chip_only")) {
      ids <- unique(d27$ENTREZID[d27$group == g])
      cat(sprintf("  GO [%s] : %d gènes", g, length(ids)))
      if (length(ids) < 10) { cat("  -> trop peu, ignoré\n"); next }
      ego <- tryCatch(enrichGO(gene = ids, universe = universe_e,
                               OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                               ont = "BP", pAdjustMethod = "BH",
                               pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                               readable = TRUE), error = function(e) NULL)
      n_term <- if (!is.null(ego)) nrow(as.data.frame(ego)) else 0
      cat(sprintf("  -> %d termes\n", n_term))
      if (n_term > 0) {
        go_res[[g]] <- as.data.frame(ego) %>% mutate(group = g)
        pdf(file.path(PLOT_DIR, sprintf("09_GO_%s_%s.pdf", FOCUS_MARK, g)),
            width = 8, height = 6)
        print(enrichplot::dotplot(ego, showCategory = 20) + ggtitle(paste(FOCUS_MARK, g)))
        dev.off()
      }
    }
    if (length(go_res))
      write_csv(bind_rows(go_res), file.path(INT_DIR, sprintf("%s_GO_by_group.csv", FOCUS_MARK)))
  } else cat("[skip] clusterProfiler absent.\n")
}

# ===========================================================================
# 10. FXN
# ===========================================================================
cat("\n=========================== 10. FXN ===========================\n")
fxn <- dat_all %>% filter(gene == "FXN") %>%
  dplyr::select(mark, type, region_class, n_regions,
                logFC_chip, FDR_chip, logFC_RNA, FDR_RNA)

if (nrow(fxn)) {
  print(as.data.frame(fxn), row.names = FALSE)
  write_csv(fxn, file.path(INT_DIR, "FXN_cross.csv"))
  cat("\nFXN n'a de pic que sur les marques narrow (H3K27ac, PolII).\n",
      "Pour les marques broad au locus FXN -> module ciblé 06/07.\n")
} else cat("FXN absent des tables. Utilise le module ciblé 06/07.\n")





# ===========================================================================
# 3b. CONCORDANCE PAR DÉCILE RNA
# ===========================================================================
binned_conc <- dat %>%
  group_by(mark) %>% mutate(bin = ntile(logFC_RNA, 10)) %>%
  group_by(mark, type, bin) %>%
  summarise(pct = round(pct_concordant(logFC_chip, logFC_RNA, type[1]), 1),
            mean_rna = mean(logFC_RNA, na.rm = TRUE),
            mean_chip = mean(logFC_chip, na.rm = TRUE),
            .groups = "drop")

print(as.data.frame(binned_conc), row.names = FALSE)
write_csv(binned_conc, file.path(INT_DIR, "concordance_by_rna_decile.csv"))

p_binned_conc <- ggplot(binned_conc, aes(bin, pct, colour = type)) +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = PAL) +
  facet_wrap(~ paste0(mark, " (", type, ")")) +
  labs(x = "Décile RNA", y = "% concordant",
       title = "Concordance par décile RNA",
       subtitle = "Est-ce que la concordance se renforce aux extrêmes ?") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

ggsave(file.path(PLOT_DIR, "03b_concordance_by_rna_decile.pdf"),
       p_binned_conc, width = 11, height = 8)



# ===========================================================================
# 6b. SCATTER AVEC FIT LINÉAIRE SIMPLIFIÉ
# ===========================================================================
p_sc2 <- ggplot(dat, aes(logFC_chip, logFC_RNA)) +
  geom_point(aes(colour = FDR_RNA < FDR_RNA_TH), alpha = 0.22, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, colour = "black") +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#4E79A7"),
                      guide = "none") +
  facet_wrap(~ paste0(mark, " (", type, ")")) +
  labs(x = "log2FC ChIP", y = "log2FC RNA",
       title = "Scatter ChIP/RNA avec fit linéaire",
       subtitle = "Lecture plus directe de la pente par marque") +
  theme_bw(base_size = 11)

ggsave(file.path(PLOT_DIR, "06b_scatter_lm.pdf"), p_sc2, width = 11, height = 8)




# ===========================================================================
# 7b. HEATMAP DES GÈNES LEADERS
# ===========================================================================
top_hm <- dat_box %>%
  group_by(mark, RNA_class) %>%
  slice_max(order_by = abs(logFC_RNA), n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = paste0(gene, " (", RNA_class, ")"))

p_hm <- ggplot(top_hm, aes(x = RNA_class, y = reorder(label, logFC_RNA), fill = logFC_chip)) +
  geom_tile() +
  facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_y") +
  scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759") +
  labs(x = NULL, y = NULL,
       title = "Top gènes par classe RNA",
       subtitle = "Lecture signature des principaux leaders") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 6))

ggsave(file.path(PLOT_DIR, "07b_top_genes_heatmap.pdf"), p_hm, width = 12, height = 8)





# ===========================================================================
# 8b. RÉSUMÉ AUC
# ===========================================================================
if (nrow(auc_tab)) {
  p_auc2 <- ggplot(auc_tab, aes(AUC, reorder(mark, AUC), fill = type)) +
    geom_col(width = 0.6) +
    geom_vline(xintercept = 0.5, linetype = 2) +
    scale_fill_manual(values = PAL) +
    labs(x = "AUC", y = NULL,
         title = "Résumé des AUC par marque",
         subtitle = "Vue compacte des performances prédictives") +
    theme_bw(base_size = 11)
  
  ggsave(file.path(PLOT_DIR, "08b_AUC_summary.pdf"), p_auc2, width = 7.5, height = 4.5)
}






# ===========================================================================
# 9b. GO EN FIGURES SÉPARÉES
# ===========================================================================
if (exists("go_res") && length(go_res)) {
  go_df <- bind_rows(go_res)
  
  for (g in unique(go_df$group)) {
    sub <- go_df %>% filter(group == g) %>%
      arrange(p.adjust) %>%
      slice_head(n = 20)
    
    p_go <- ggplot(sub, aes(x = GeneRatio, y = reorder(Description, GeneRatio),
                            size = Count, colour = p.adjust)) +
      geom_point() +
      scale_colour_gradient(low = "#E15759", high = "#4E79A7", trans = "reverse") +
      labs(x = "Gene ratio", y = NULL,
           title = paste("GO", FOCUS_MARK, g),
           subtitle = "Top 20 terms, one plot per group") +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_text(size = 7))
    
    ggsave(file.path(PLOT_DIR, sprintf("09b_GO_%s_%s_summary.pdf", FOCUS_MARK, g)),
           p_go, width = 9, height = 6)
  }
  
  write_csv(go_df, file.path(INT_DIR, sprintf("%s_GO_by_group.csv", FOCUS_MARK)))
}




# ===========================================================================
# 10b. RÉSUMÉ FXN
# ===========================================================================
fxn_zoom <- dat_all %>% filter(gene == "FXN")

if (nrow(fxn_zoom)) {
  p_fxn <- ggplot(fxn_zoom, aes(mark, logFC_chip, fill = type)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("RNA %.2f", logFC_RNA)), vjust = -0.2, size = 3) +
    scale_fill_manual(values = PAL) +
    labs(x = NULL, y = "FXN log2FC ChIP",
         title = "Résumé locus FXN",
         subtitle = "Synthèse single-gene across marks") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  
  ggsave(file.path(PLOT_DIR, "10b_FXN_summary.pdf"), p_fxn, width = 8, height = 4.5)
}



# ===========================================================================
# 11. MATRICE DE CORRÉLATION
# ===========================================================================
cm <- dat_all %>%
  dplyr::select(logFC_chip, logFC_RNA, baseMean, n_regions, mean_signal) %>%
  mutate(across(everything(), as.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "r")

p_cm <- ggplot(cm, aes(var1, var2, fill = r)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759", limits = c(-1, 1)) +
  labs(x = NULL, y = NULL,
       title = "Matrice de corrélation étendue",
       subtitle = "Résumé global des variables d'intégration") +
  theme_bw(base_size = 11)

ggsave(file.path(PLOT_DIR, "11b_correlation_matrix.pdf"), p_cm, width = 8, height = 6)

# ===========================================================================
# 12. BOXPLOT PAR CLASSE DE RÉGION
# ===========================================================================
p_reg_box <- ggplot(dat_all %>% filter(region_class %in% c("Promoter", "Genic", "Distal")),
                    aes(region_class, logFC_chip, fill = region_class)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~ paste0(mark, " (", type, ")")) +
  labs(x = NULL, y = "log2FC ChIP",
       title = "Effet ChIP par classe de région",
       subtitle = "Promoter vs genic vs distal") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

ggsave(file.path(PLOT_DIR, "12b_chip_by_region_boxplot.pdf"), p_reg_box, width = 11, height = 8)












# ==========================================================================
# FA-CENTRIC CHUNKS — version corrigée sans GeneSymbol
# ========================================================================== 

# ==========================================================================
# 9c. FA PATHWAY SCORES — résumé par modules biologiques
# ==========================================================================
fa_modules <- list(
  Mitochondria = c("FXN","NDUFS1","NDUFS2","NDUFA1","NDUFA2","NDUFB8","SDHA","SDHB","UQCRC1","UQCRC2","COX4I1","COX5A","ATP5F1A"),
  OxStress = c("SOD1","SOD2","GPX1","GPX4","CAT","PRDX1","PRDX3","TXN","TXNRD1","HMOX1","NQO1"),
  Neuronal = c("NEFL","NEFM","MAP2","TUBB3","RBFOX3","SYN1","SLC17A7","GAD1","GAD2","STMN2"),
  Cardiac = c("TNNT2","TTN","MYH6","MYH7","ACTC1","RYR2","CACNA1C","PLN","ATP2A2"),
  Adhesion = c("ITGA6","ITGB1","CDH1","CDH2","ICAM1","VCAM1","FN1","COL1A1","COL1A2")
)

fa_df <- dat_all %>%
  mutate(symbol = if ("GeneSymbol" %in% names(.)) coalesce(GeneSymbol, gene) else gene) %>%
  filter(!is.na(symbol)) %>%
  mutate(module = case_when(
    symbol %in% fa_modules$Mitochondria ~ "Mitochondria",
    symbol %in% fa_modules$OxStress     ~ "OxStress",
    symbol %in% fa_modules$Neuronal     ~ "Neuronal",
    symbol %in% fa_modules$Cardiac      ~ "Cardiac",
    symbol %in% fa_modules$Adhesion     ~ "Adhesion",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(module))

fa_score <- fa_df %>%
  group_by(mark, type, module) %>%
  summarise(
    n = n(),
    mean_rna = mean(logFC_RNA, na.rm = TRUE),
    mean_chip = mean(logFC_chip, na.rm = TRUE),
    med_rna = median(logFC_RNA, na.rm = TRUE),
    med_chip = median(logFC_chip, na.rm = TRUE),
    .groups = "drop"
  )

print(as.data.frame(fa_score), row.names = FALSE)
write_csv(fa_score, file.path(INT_DIR, "FA_module_scores.csv"))

p_fa_score <- ggplot(fa_score, aes(module, mean_chip, fill = type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  scale_fill_manual(values = PAL) +
  facet_wrap(~ paste0(mark, " (", type, ")")) +
  labs(x = NULL, y = "Mean log2FC ChIP",
       title = "FA module scores",
       subtitle = "Mean chromatin effect by biological module") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(PLOT_DIR, "09c_FA_module_scores.pdf"), p_fa_score, width = 11, height = 8)

# ==========================================================================
# 9d. HEATMAP FA — gènes modules + effets RNA/ChIP
# ==========================================================================
fa_heat <- fa_df %>%
  mutate(label = symbol) %>%
  select(symbol, module, mark, type, logFC_RNA, logFC_chip, FDR_RNA, FDR_chip) %>%
  distinct()

p_fa_heat <- ggplot(fa_heat, aes(x = mark, y = reorder(symbol, logFC_RNA), fill = logFC_chip)) +
  geom_tile() +
  facet_grid(module ~ type, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759") +
  labs(x = NULL, y = NULL,
       title = "FA gene heatmap",
       subtitle = "Selected FA-relevant genes across marks and module classes") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 6))

ggsave(file.path(PLOT_DIR, "09d_FA_gene_heatmap.pdf"), p_fa_heat, width = 12, height = 10)

# ==========================================================================
# 9e. SCATTER FA — RNA vs ChIP sur gènes sélectionnés
# ==========================================================================
module_genes <- unique(unlist(fa_modules))
fa_scatter_df <- dat_all %>%
  mutate(symbol = if ("GeneSymbol" %in% names(.)) coalesce(GeneSymbol, gene) else gene) %>%
  filter(symbol %in% module_genes) %>%
  mutate(module = case_when(
    symbol %in% fa_modules$Mitochondria ~ "Mitochondria",
    symbol %in% fa_modules$OxStress     ~ "OxStress",
    symbol %in% fa_modules$Neuronal     ~ "Neuronal",
    symbol %in% fa_modules$Cardiac      ~ "Cardiac",
    symbol %in% fa_modules$Adhesion     ~ "Adhesion",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(module))

p_fa_scatter <- ggplot(fa_scatter_df, aes(logFC_chip, logFC_RNA)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(colour = module), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8, colour = "black", na.rm = TRUE) +
  ggrepel::geom_text_repel(
    data = fa_scatter_df %>% filter(symbol == "FXN"),
    aes(label = symbol),
    size = 3, max.overlaps = 20
  ) +
  facet_wrap(~ paste0(mark, " (", type, ")")) +
  scale_colour_brewer(palette = "Set2") +
  labs(x = "log2FC ChIP", y = "log2FC RNA",
       title = "FA-focused scatter",
       subtitle = "Selected pathway genes only") +
  theme_bw(base_size = 11)

ggsave(file.path(PLOT_DIR, "09e_FA_scatter.pdf"), p_fa_scatter, width = 11, height = 8)

# ==========================================================================
# 10c. FXN AND NEIGHBOURS
# ==========================================================================
fxn_set <- c("FXN","NFS1","ISCU","LIAS","FDX2","NFU1","BOLA1","BOLA3","LIPT1","LIPT2","GLRX5")

fxn_nei <- dat_all %>%
  mutate(symbol = if ("GeneSymbol" %in% names(.)) coalesce(GeneSymbol, gene) else gene) %>%
  filter(symbol %in% fxn_set) %>%
  select(symbol, mark, type, logFC_RNA, logFC_chip) %>%
  distinct()

if (nrow(fxn_nei)) {
  p_fxn_nei <- ggplot(fxn_nei, aes(symbol, logFC_chip, fill = type)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    scale_fill_manual(values = PAL) +
    facet_wrap(~ paste0(mark, " (", type, ")")) +
    labs(x = NULL, y = "log2FC ChIP",
         title = "FXN iron-sulfur module",
         subtitle = "Core FA / mitochondrial genes") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave(file.path(PLOT_DIR, "10c_FXN_module.pdf"), p_fxn_nei, width = 12, height = 7)
}

# ==========================================================================
# 9f. GROUPES FA x CONCORDANCE
# ==========================================================================
fa_group <- dat_all %>%
  mutate(symbol = if ("GeneSymbol" %in% names(.)) coalesce(GeneSymbol, gene) else gene,
         fc_rna_dir = case_when(FDR_RNA < FDR_RNA_TH & logFC_RNA > 0 ~ "RNA_up",
                                FDR_RNA < FDR_RNA_TH & logFC_RNA < 0 ~ "RNA_down",
                                TRUE ~ "RNA_ns"),
         fc_chip_dir = case_when(logFC_chip > 0 ~ "ChIP_up", logFC_chip < 0 ~ "ChIP_down", TRUE ~ "ChIP_0")) %>%
  mutate(concord = case_when(
    (logFC_chip > 0 & logFC_RNA > 0) ~ "Concordant_up",
    (logFC_chip < 0 & logFC_RNA < 0) ~ "Concordant_down",
    TRUE ~ "Discordant"
  )) %>%
  filter(symbol %in% module_genes)

fa_group_sum <- fa_group %>%
  count(concord, name = "n")

write_csv(fa_group_sum, file.path(INT_DIR, "FA_concordance_groups.csv"))

p_fa_group <- ggplot(fa_group_sum, aes(concord, n, fill = concord)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(Concordant_up = "#E15759", Concordant_down = "#4E79A7", Discordant = "grey70")) +
  labs(x = NULL, y = "n",
       title = "FA genes by concordance class",
       subtitle = "How many FA-relevant genes are concordant vs discordant ?") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(PLOT_DIR, "09f_FA_concordance_groups.pdf"), p_fa_group, width = 8, height = 5)

# ==========================================================================
# 9g. GO — ONE GROUP / ONE PLOT
# ==========================================================================
if (exists("go_res") && length(go_res)) {
  go_df <- bind_rows(go_res)
  
  for (g in unique(go_df$group)) {
    sub <- go_df %>%
      filter(group == g) %>%
      arrange(p.adjust) %>%
      slice_head(n = 15)
    
    p_go <- ggplot(sub, aes(x = GeneRatio, y = reorder(Description, GeneRatio),
                            size = Count, colour = p.adjust)) +
      geom_point() +
      scale_colour_gradient(low = "#E15759", high = "#4E79A7", trans = "reverse") +
      labs(x = "Gene ratio", y = NULL,
           title = paste("GO", FOCUS_MARK, g),
           subtitle = "Top 15 terms, separate plots per group") +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_text(size = 7))
    
    ggsave(file.path(PLOT_DIR, sprintf("09g_GO_%s_%s.pdf", FOCUS_MARK, g)),
           p_go, width = 9, height = 6)
  }
}











# ===========================================================================
# RÉCAP
# ===========================================================================
summary_tab <- conc %>%
  dplyr::select(mark, type, n, pct_observed = pct, OR, fisher_p) %>%
  left_join(perm %>% dplyr::select(mark, null_mean, excess, z, p_perm), by = "mark")
if (nrow(auc_tab))
  summary_tab <- summary_tab %>%
  left_join(auc_tab %>% dplyr::select(mark, AUC_signed = AUC), by = "mark")

cat("\n=========================== RÉCAP ===========================\n")
print(as.data.frame(summary_tab %>% mutate(fisher_p = signif(fisher_p, 3))),
      row.names = FALSE)
write_csv(summary_tab, file.path(INT_DIR, "SUMMARY.csv"))

cat("\n==> Figures : ", PLOT_DIR, "\n\n")





