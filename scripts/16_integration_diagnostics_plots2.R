#!/usr/bin/env Rscript
# =============================================================================
# 16_integration_analysis.R — Intégration RNA-seq x ChIP-seq (FRDA)   [v3]
# -----------------------------------------------------------------------------
# PRÉREQUIS : 04b_recompute_gene_level.R
#
# CORRECTIFS v3 (par rapport à ta version)
#  [BUG] concordance FA ignorait le type de marque -> "même signe = concordant"
#        est FAUX pour H3K9me3/H3K27me3. Corrigé via concord_class().
#  [BUG] les sections FA utilisaient dat_all (1 ligne par gène x classe de
#        région) -> gènes comptés jusqu'à 4x, moyennes pondérées par le nombre
#        de classes, geom_col empilant 4 barres. Basculé sur `dat`.
#  [BUG] matrice de corrélation : logFC_chip des 4 marques mélangés dans une
#        seule colonne -> r sans signification. Corrigé par pivot large.
#  [BUG] deux dplyr::select non qualifiés (§9d, §10c) -> plantage S4.
#  [BUG] GeneRatio est une chaîne "5/100" -> axe X illisible. Parsé.
#  [BIO] modules Neuronal/Cardiac non exprimés en fibroblaste -> remplacés par
#        des modules FRDA (Fe-S, fer, OXPHOS, stress ox, biogenèse mito) +
#        un module de CONTRÔLE NÉGATIF explicite.
#  [NEW] CHIP_AGG : bascule mean_logFC / max_abs pour tester la robustesse.
#  [NEW] §0d : test de robustesse de l'agrégation.
#  - §9b et §9g étaient redondants -> fusionnés en une fonction unique.
#
# Usage : export PROJECT_ROOT=$(pwd) && Rscript 16_integration_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(tibble)
  library(ggplot2); library(ggrepel); library(stringr)
  library(org.Hs.eg.db); library(AnnotationDbi)
})

set.seed(42)

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
ANN_DIR  <- file.path(PROJECT_ROOT, "results", "annotation")
INT_DIR  <- file.path(PROJECT_ROOT, "results", "integration")
PLOT_DIR <- file.path(INT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

RNA_FILE <- "/home/apuhart/Downloads/Bulk-RNAseq/results_test/fibro/control_vs_FA/fibro_control_vs_FA_DESeq2_results.csv"

MARK_TYPE <- c(H3K9me3 = "repressive", H3K27me3 = "repressive",
               H3K27ac = "active",     PolII    = "active")

# AGRÉGATION DES RÉGIONS EN GÈNE — paramètre critique.
#   "maxabs" = colonne logFC  -> région au |fold| le plus fort. Sur des folds
#              NON shrunkés, cela sélectionne activement la région la MOINS
#              couverte (la plus bruitée) de chaque gène.
#   "mean"   = colonne mean_logFC -> moyenne des régions. Plus robuste.
# Lance les deux et compare : si la concordance tient, le résultat est robuste.
CHIP_AGG <- "mean"          # "mean" | "maxabs"

N_PERM      <- 1000
FDR_RNA_TH  <- 0.05
LFC_RNA_TH  <- 1
LFC_RNA_MAX <- 15

PAL <- c(repressive = "#4E79A7", active = "#E15759")

save_plot <- function(p, filename, width = 11, height = 8) {
  ggsave(file.path(PLOT_DIR, filename), p, width = width, height = height)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Helpers de concordance — TYPE-AWARE (le bug n°1 de la version précédente)
# ---------------------------------------------------------------------------
# répressive : concordant = signes OPPOSÉS (plus de marque -> moins d'ARN)
# active     : concordant = MÊME signe
concord_class <- function(lfc_chip, lfc_rna, type) {
  prod <- lfc_chip * lfc_rna
  ifelse(is.na(prod) | prod == 0, NA_character_,
         ifelse((type == "active" & prod > 0) | (type == "repressive" & prod < 0),
                "Concordant", "Discordant"))
}

pct_concordant <- function(lfc_chip, lfc_rna, type) {
  cc <- concord_class(lfc_chip, lfc_rna, type)
  cc <- cc[!is.na(cc)]
  if (!length(cc)) return(NA_real_)
  100 * mean(cc == "Concordant")
}

# ---------------------------------------------------------------------------
# 0a. RNA-seq
# ---------------------------------------------------------------------------
cat("\n=========================== CHARGEMENT ===========================\n")
cat(sprintf("Agrégation ChIP : %s\n", CHIP_AGG))

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

chip_raw <- lapply(files, function(f) {
  tag <- sub("_gene_level_fixed$", "", tools::file_path_sans_ext(basename(f)))
  read_csv(f, show_col_types = FALSE) %>%
    mutate(mark = sub("__.*$", "", tag), ENTREZID = as.character(GeneID))
}) %>% bind_rows() %>%
  distinct(mark, ENTREZID, region_class, .keep_all = TRUE)

if (!"mean_logFC" %in% names(chip_raw)) {
  warning("Colonne mean_logFC absente -> repli sur logFC (max_abs).")
  chip_raw$mean_logFC <- chip_raw$logFC
}

chip <- chip_raw %>%
  mutate(logFC_chip = if (CHIP_AGG == "mean") mean_logFC else logFC)

build_dat <- function(chip_df, agg) {
  chip_df %>%
    mutate(logFC_chip = if (agg == "mean") mean_logFC else logFC) %>%
    inner_join(rna, by = "ENTREZID", relationship = "many-to-one") %>%
    transmute(mark, type = unname(MARK_TYPE[mark]), region_class,
              gene = coalesce(GeneSymbol, SYMBOL),
              ENTREZID, n_regions, mean_signal,
              logFC_chip, FDR_chip = FDR,
              baseMean, logFC_RNA, PValue_RNA, FDR_RNA) %>%
    filter(!is.na(logFC_chip), !is.na(logFC_RNA))
}

dat_all <- build_dat(chip_raw, CHIP_AGG)

# `dat` = UNE ligne par gène x marque. C'est la vue à utiliser partout sauf
# dans les analyses explicitement stratifiées par classe de région (§4, §12).
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

# ---------------------------------------------------------------------------
# 0d. ROBUSTESSE DE L'AGRÉGATION  <<< NOUVEAU
# ---------------------------------------------------------------------------
# "maxabs" sélectionne la région au fold le plus extrême. Sur des folds non
# shrunkés, c'est souvent la région la moins couverte (la plus bruitée).
# Si la concordance change beaucoup entre "mean" et "maxabs", le résultat
# dépend de la région choisie -> fragile. Si elle tient -> robuste.
cat("\n=================== 0d. ROBUSTESSE DE L'AGRÉGATION ===================\n")
robust <- lapply(c("mean", "maxabs"), function(agg) {
  d <- build_dat(chip_raw, agg) %>% filter(region_class == "All")
  d %>% group_by(mark, type) %>%
    summarise(pct = round(pct_concordant(logFC_chip, logFC_RNA, type[1]), 1),
              .groups = "drop") %>%
    mutate(agg = agg)
}) %>% bind_rows() %>%
  pivot_wider(names_from = agg, values_from = pct) %>%
  mutate(delta = round(mean - maxabs, 1))
print(as.data.frame(robust), row.names = FALSE)
cat("\ndelta = concordance(mean) - concordance(maxabs).\n")
cat("|delta| < 2 points -> robuste. |delta| > 5 -> le résultat dépend de\n")
cat("la région choisie : à signaler, et préférer 'mean'.\n")
write_csv(robust, file.path(INT_DIR, "aggregation_robustness.csv"))

# ===========================================================================
# 1. CONCORDANCE GLOBALE
# ===========================================================================
# Le % brut dépend des marges ; l'ODDS RATIO n'en dépend pas.
# Sens attendu : OR > 1 si ACTIVE, OR < 1 si RÉPRESSIVE.
cat("\n=================== 1. CONCORDANCE GLOBALE ===================\n")

conc <- lapply(split(dat, dat$mark), function(d) {
  d <- d %>% filter(logFC_chip != 0, logFC_RNA != 0)
  tab <- table(chip = sign(d$logFC_chip), rna = sign(d$logFC_RNA))
  if (!all(dim(tab) == c(2, 2)))
    return(tibble(mark = d$mark[1], type = d$type[1], n = nrow(d),
                  OR = NA_real_, fisher_p = NA_real_, pct = NA_real_))
  ft <- fisher.test(tab)
  tibble(mark = d$mark[1], type = d$type[1], n = sum(tab),
         OR = round(unname(ft$estimate), 3), fisher_p = ft$p.value,
         pct = round(pct_concordant(d$logFC_chip, d$logFC_RNA, d$type[1]), 1))
}) %>% bind_rows() %>% arrange(desc(pct))

print(as.data.frame(conc %>% mutate(fisher_p = signif(fisher_p, 3))), row.names = FALSE)

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
         observed = round(obs, 2), null_mean = round(mean(nv), 2),
         excess = round(obs - mean(nv), 2), null_sd = round(sd(nv), 3),
         z = round((obs - mean(nv)) / sd(nv), 1),
         p_perm = (1 + sum(abs(nv - mean(nv)) >= abs(obs - mean(nv)))) / (N_PERM + 1))
}) %>% bind_rows() %>% arrange(desc(excess))

print(as.data.frame(perm), row.names = FALSE)
write_csv(perm, file.path(INT_DIR, "permutation_null.csv"))
cat("\nexcess = TAILLE de l'effet | z = SIGNIFICATIVITÉ.\n")
cat("Le null n'est PAS à 50% : le comparer à 50% donne des conclusions fausses.\n")

perm_plot <- lapply(split(dat, dat$mark), function(d) {
  d <- d %>% filter(logFC_chip != 0, logFC_RNA != 0)
  ty <- d$type[1]
  tibble(mark = d$mark[1], type = ty,
         null = replicate(300, pct_concordant(sample(d$logFC_chip), d$logFC_RNA, ty)))
}) %>% bind_rows() %>% left_join(perm %>% dplyr::select(mark, observed), by = "mark")

save_plot(
  ggplot(perm_plot, aes(null)) +
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
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "02_permutation_null.pdf")

# ===========================================================================
# 3. STRATIFICATION PAR ABONDANCE  <<< LE CONTRÔLE DÉCISIF
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

save_plot(
  ggplot(strat, aes(factor(q_signal), pct, fill = type)) +
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
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "03_concordance_by_abundance.pdf", 10, 7)

# ---- 3b. Concordance par décile RNA ---------------------------------------
binned_conc <- dat %>%
  group_by(mark) %>% mutate(bin = ntile(logFC_RNA, 10)) %>%
  group_by(mark, type, bin) %>%
  summarise(n = n(),
            pct = round(pct_concordant(logFC_chip, logFC_RNA, type[1]), 1),
            mean_rna = round(mean(logFC_RNA, na.rm = TRUE), 3),
            mean_chip = round(mean(logFC_chip, na.rm = TRUE), 3),
            .groups = "drop")
write_csv(binned_conc, file.path(INT_DIR, "concordance_by_rna_decile.csv"))

save_plot(
  ggplot(binned_conc, aes(bin, pct, colour = type)) +
    geom_hline(yintercept = 50, linetype = 2, colour = "grey60") +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    scale_x_continuous(breaks = 1:10) +
    scale_colour_manual(values = PAL) +
    facet_wrap(~ paste0(mark, " (", type, ")")) +
    labs(x = "Décile RNA (1 = plus réprimé → 10 = plus induit)", y = "% concordant",
         title = "Concordance par décile RNA",
         subtitle = "La concordance se renforce-t-elle aux extrêmes ?") +
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "03b_concordance_by_rna_decile.pdf")

# ===========================================================================
# 4. PROMOTEUR vs DISTAL   (utilise dat_all : c'est le but)
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
  save_plot(
    ggplot(by_reg, aes(region_class, pct, fill = type)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 50, linetype = 2, colour = "grey30") +
      geom_text(aes(label = paste0(pct, "%\nn=", n)), vjust = -0.3, size = 2.6) +
      scale_fill_manual(values = PAL) +
      facet_wrap(~ paste0(mark, " (", type, ")")) +
      coord_cartesian(ylim = c(0, 85)) +
      labs(x = NULL, y = "% concordant",
           title = "Concordance par classe de région",
           subtitle = "Promoteur = régulation directe | Distal = paysage enhancer") +
      theme_bw(base_size = 11) + theme(legend.position = "none"),
    "04_concordance_by_region.pdf", 10, 7)
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

save_plot(
  ggplot(binned, aes(bin, mean_chip)) +
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
         subtitle = "Répressive : pente descendante | Active : pente montante") +
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "05_binned_trend.pdf")

# ===========================================================================
# 6. SCATTER
# ===========================================================================
lab <- dat %>% group_by(mark) %>%
  filter(gene == "FXN" | rank(FDR_RNA, ties.method = "first") <= 8) %>% ungroup()

save_plot(
  ggplot(dat, aes(logFC_chip, logFC_RNA)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(aes(colour = FDR_RNA < FDR_RNA_TH), alpha = 0.35, size = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "#E15759", linewidth = 0.7) +
    geom_text_repel(data = lab, aes(label = gene), size = 2.4,
                    max.overlaps = 20, min.segment.length = 0) +
    geom_point(data = filter(dat, gene == "FXN"),
               colour = "black", size = 3, shape = 21, stroke = 1.1, fill = NA) +
    scale_colour_manual(values = c(`FALSE` = "grey78", `TRUE` = "#4E79A7"),
                        name = paste0("RNA FDR < ", FDR_RNA_TH)) +
    facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_x") +
    labs(x = "log2FC ChIP (FRDA vs Control)", y = "log2FC RNA (FRDA vs Control)",
         title = "Chromatine vs expression", subtitle = "Cercle noir = FXN") +
    theme_bw(base_size = 11),
  "06_scatter.pdf")

# ===========================================================================
# 7. BOXPLOT + Kruskal-Wallis
# ===========================================================================
dat_box <- dat %>%
  mutate(RNA_class = factor(case_when(
    FDR_RNA < FDR_RNA_TH & logFC_RNA > 0 ~ "Up",
    FDR_RNA < FDR_RNA_TH & logFC_RNA < 0 ~ "Down",
    TRUE ~ "Stable"), levels = c("Down", "Stable", "Up")))

n_lab <- dat_box %>% count(mark, type, RNA_class)

save_plot(
  ggplot(dat_box, aes(RNA_class, logFC_chip, fill = RNA_class)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_boxplot(outlier.size = 0.2, notch = TRUE) +
    geom_text(data = n_lab, aes(x = RNA_class, y = -Inf, label = paste0("n=", n)),
              vjust = -0.6, size = 2.5, inherit.aes = FALSE) +
    facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_y") +
    scale_fill_manual(values = c(Down = "#4E79A7", Stable = "grey85", Up = "#E15759")) +
    labs(x = "Catégorie RNA", y = "log2FC ChIP",
         title = "Changement chromatinien selon l'état transcriptionnel",
         subtitle = "Encoches disjointes = médianes significativement différentes") +
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "07_boxplot_by_RNA_class.pdf", 10, 7)

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
    save_plot(
      ggplot(auc_tab, aes(reorder(mark, AUC), AUC, fill = type)) +
        geom_col(width = 0.6) +
        geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15) +
        geom_hline(yintercept = 0.5, linetype = 2) +
        geom_text(aes(label = paste0(AUC, "\nn=", n)), hjust = -0.2, size = 3) +
        scale_fill_manual(values = PAL) +
        coord_flip(ylim = c(0.3, 1)) +
        labs(x = NULL, y = "ROC AUC (± IC95)",
             title = "La chromatine prédit-elle le SENS du changement RNA ?") +
        theme_bw(base_size = 11),
      "08_AUC_signed.pdf", 7, 4.5)
  }
} else cat("\n[skip] package pROC absent.\n")

# ===========================================================================
# 9. PARTITION DIRECT / INDIRECT + GO
# ===========================================================================
cat("\n=================== 9. DIRECT / INDIRECT ===================\n")

FOCUS_MARK <- "H3K27ac"
d27 <- dat %>% filter(mark == FOCUS_MARK)
go_res <- list()

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
      if (n_term > 0) go_res[[g]] <- as.data.frame(ego) %>% mutate(group = g)
    }
  } else cat("[skip] clusterProfiler absent.\n")
}

# ---- 9b. GO : une figure par groupe (fusion des anciens 9b et 9g) ---------
# enrichGO renvoie GeneRatio sous forme de chaîne "5/100" -> il FAUT le parser,
# sinon ggplot en fait un facteur et l'axe X est illisible.
if (length(go_res)) {
  go_df <- bind_rows(go_res) %>%
    mutate(GeneRatio_num = vapply(strsplit(GeneRatio, "/"),
                                  function(x) as.numeric(x[1]) / as.numeric(x[2]),
                                  numeric(1)))
  write_csv(go_df, file.path(INT_DIR, sprintf("%s_GO_by_group.csv", FOCUS_MARK)))
  
  for (g in unique(go_df$group)) {
    sub <- go_df %>% filter(group == g) %>% arrange(p.adjust) %>% slice_head(n = 20)
    save_plot(
      ggplot(sub, aes(GeneRatio_num, reorder(Description, GeneRatio_num),
                      size = Count, colour = p.adjust)) +
        geom_point() +
        scale_colour_gradient(low = "#E15759", high = "#4E79A7") +
        labs(x = "Gene ratio", y = NULL,
             title = paste("GO:BP —", FOCUS_MARK, "—", g),
             subtitle = "Top 20 termes") +
        theme_bw(base_size = 11) + theme(axis.text.y = element_text(size = 7)),
      sprintf("09b_GO_%s_%s.pdf", FOCUS_MARK, g), 9, 6)
  }
}

# ===========================================================================
# 9c-9f. MODULES FRDA
# ===========================================================================
# Modules révisés : les anciens "Neuronal" et "Cardiac" ne sont PAS exprimés
# dans un fibroblaste -> ils ne mesuraient que du bruit de comptages nuls.
# Ils sont conservés sous le nom NegCtrl_notExpressed, et servent désormais
# de CONTRÔLE NÉGATIF explicite : on n'attend AUCUN signal dessus. S'il y en a,
# c'est un artefact ailleurs dans la chaîne.
fa_modules <- list(
  FeS_assembly = c("FXN","NFS1","ISCU","ISCA1","ISCA2","FDX1","FDX2","NUBPL",
                   "NFU1","BOLA1","BOLA3","LYRM4","GLRX5","LIAS","LIPT1","LIPT2",
                   "HSPA9","HSCB"),
  Iron_homeostasis = c("TFRC","SLC40A1","FTH1","FTL","SLC11A2","ACO1","ACO2","IREB2"),
  OXPHOS = c("NDUFS1","NDUFS2","NDUFA9","NDUFB8","SDHA","SDHB","UQCRC1","UQCRC2",
             "COX4I1","COX5A","ATP5F1A"),
  Oxidative_stress = c("NFE2L2","HMOX1","SOD1","SOD2","GPX1","GPX4","NQO1","GCLC",
                       "CAT","PRDX1","PRDX3","TXN","TXNRD1"),
  Mito_biogenesis = c("PPARGC1A","TFAM","NRF1","ESRRA","PPRC1"),
  NegCtrl_notExpressed = c("NEFL","NEFM","MAP2","TUBB3","RBFOX3","SYN1","GAD1",
                           "TNNT2","MYH6","MYH7","RYR2")
)
fa_lookup <- tibble(gene = unlist(fa_modules, use.names = FALSE),
                    module = rep(names(fa_modules), lengths(fa_modules)))

cat("\n=================== 9c. MODULES FRDA ===================\n")

# `dat` et non `dat_all` : sinon chaque gène est compté jusqu'à 4x (une fois
# par classe de région) et les moyennes sont pondérées par le nombre de classes.
fa_df <- dat %>% inner_join(fa_lookup, by = "gene")

if (nrow(fa_df)) {
  fa_score <- fa_df %>%
    group_by(mark, type, module) %>%
    summarise(n = n(),
              mean_rna  = round(mean(logFC_RNA, na.rm = TRUE), 3),
              mean_chip = round(mean(logFC_chip, na.rm = TRUE), 3),
              pct_conc  = round(pct_concordant(logFC_chip, logFC_RNA, type[1]), 1),
              .groups = "drop")
  print(as.data.frame(fa_score), row.names = FALSE)
  write_csv(fa_score, file.path(INT_DIR, "FA_module_scores.csv"))
  cat("\nNegCtrl_notExpressed : contrôle négatif. Aucun signal attendu.\n")
  
  save_plot(
    ggplot(fa_score, aes(module, mean_chip, fill = type)) +
      geom_col(width = 0.65) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
      geom_text(aes(label = paste0("n=", n)), vjust = -0.3, size = 2.3) +
      scale_fill_manual(values = PAL) +
      facet_wrap(~ paste0(mark, " (", type, ")")) +
      labs(x = NULL, y = "log2FC ChIP moyen",
           title = "Modules FRDA — effet chromatinien moyen",
           subtitle = "NegCtrl = contrôle négatif (gènes non exprimés en fibroblaste)") +
      theme_bw(base_size = 11) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 40, hjust = 1, size = 7)),
    "09c_FA_module_scores.pdf")
  
  # ---- 9d. Heatmap des gènes de modules ----------------------------------
  save_plot(
    ggplot(fa_df, aes(mark, reorder(gene, logFC_RNA), fill = logFC_chip)) +
      geom_tile(colour = "white") +
      facet_grid(module ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759") +
      labs(x = NULL, y = NULL, fill = "log2FC\nChIP",
           title = "Gènes FRDA à travers les marques",
           subtitle = "Ordonnés par log2FC RNA") +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_text(size = 6),
            strip.text.y = element_text(size = 6, angle = 0)),
    "09d_FA_gene_heatmap.pdf", 9, 12)
  
  # ---- 9e. Scatter FRDA ---------------------------------------------------
  save_plot(
    ggplot(fa_df, aes(logFC_chip, logFC_RNA)) +
      geom_hline(yintercept = 0, colour = "grey85") +
      geom_vline(xintercept = 0, colour = "grey85") +
      geom_point(aes(colour = module), alpha = 0.75, size = 2) +
      geom_text_repel(data = fa_df %>% filter(gene == "FXN"), aes(label = gene),
                      size = 3, max.overlaps = 20) +
      facet_wrap(~ paste0(mark, " (", type, ")")) +
      scale_colour_brewer(palette = "Set2") +
      labs(x = "log2FC ChIP", y = "log2FC RNA",
           title = "Scatter restreint aux modules FRDA") +
      theme_bw(base_size = 11),
    "09e_FA_scatter.pdf")
  
  # ---- 9f. Concordance des gènes FRDA  [BUG CORRIGÉ] ---------------------
  # L'ancienne version étiquetait "même signe = concordant" pour TOUTES les
  # marques : faux pour H3K9me3/H3K27me3, où concordant = signes opposés.
  fa_group <- fa_df %>%
    mutate(concord = concord_class(logFC_chip, logFC_RNA, type)) %>%
    filter(!is.na(concord))
  
  fa_group_sum <- fa_group %>% count(mark, type, module, concord)
  write_csv(fa_group_sum, file.path(INT_DIR, "FA_concordance_groups.csv"))
  
  save_plot(
    ggplot(fa_group_sum, aes(module, n, fill = concord)) +
      geom_col(position = "fill", width = 0.7) +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "grey30") +
      scale_fill_manual(values = c(Concordant = "#59A14F", Discordant = "grey75")) +
      scale_y_continuous(labels = scales::percent) +
      facet_wrap(~ paste0(mark, " (", type, ")")) +
      labs(x = NULL, y = "Proportion", fill = NULL,
           title = "Concordance des gènes FRDA, par module",
           subtitle = "Concordant = sens attendu selon le TYPE de marque") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7)),
    "09f_FA_concordance.pdf")
} else cat("Aucun gène de module FRDA trouvé dans les données.\n")

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
  
  # `dat` : une barre par marque. Avec dat_all, geom_col empilerait les 4
  # classes de région sans le dire.
  fxn_one <- dat %>% filter(gene == "FXN")
  if (nrow(fxn_one))
    save_plot(
      ggplot(fxn_one, aes(mark, logFC_chip, fill = type)) +
        geom_col(width = 0.65) +
        geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
        geom_text(aes(label = sprintf("ChIP %.3f\nRNA %.2f", logFC_chip, logFC_RNA)),
                  vjust = -0.2, size = 3) +
        scale_fill_manual(values = PAL) +
        labs(x = NULL, y = "log2FC ChIP",
             title = "Locus FXN à travers les marques",
             subtitle = sprintf("RNA log2FC = %.2f (FDR = %.3g)",
                                fxn_one$logFC_RNA[1], fxn_one$FDR_RNA[1])) +
        theme_bw(base_size = 11) + theme(legend.position = "none"),
      "10b_FXN_summary.pdf", 8, 4.5)
  
  cat("\nFXN n'a de pic que sur les marques narrow. Pour les broad au locus FXN\n",
      "-> module ciblé 06/07 (l'assignation pic->gène ne convient pas aux domaines).\n")
} else cat("FXN absent des tables. Utilise le module ciblé 06/07.\n")

# ---- 10c. Module Fe-S ------------------------------------------------------
fes <- dat %>% filter(gene %in% fa_modules$FeS_assembly) %>%
  dplyr::select(gene, mark, type, logFC_RNA, FDR_RNA, logFC_chip)

if (nrow(fes)) {
  write_csv(fes, file.path(INT_DIR, "FeS_module_cross.csv"))
  save_plot(
    ggplot(fes, aes(reorder(gene, logFC_chip), logFC_chip, fill = type)) +
      geom_col(width = 0.65) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
      scale_fill_manual(values = PAL) +
      facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_x") +
      labs(x = NULL, y = "log2FC ChIP",
           title = "Module d'assemblage Fe-S",
           subtitle = "Cœur de la cascade frataxine") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "none"),
    "10c_FeS_module.pdf", 12, 7)
}

# ===========================================================================
# 11. CORRÉLATION ENTRE MARQUES  [BUG CORRIGÉ]
# ===========================================================================
# L'ancienne version corrélait le logFC_chip des 4 marques MÉLANGÉES dans une
# seule colonne, sur des lignes dupliquées 4x -> le r n'avait aucun sens.
# Ici : pivot large, une colonne par marque, une ligne par gène.
cat("\n=================== 11. CORRÉLATION ENTRE MARQUES ===================\n")

wide <- dat %>%
  dplyr::select(ENTREZID, mark, logFC_chip) %>%
  pivot_wider(names_from = mark, values_from = logFC_chip) %>%
  left_join(dat %>% distinct(ENTREZID, logFC_RNA), by = "ENTREZID")

cor_mat <- wide %>% dplyr::select(-ENTREZID) %>%
  cor(use = "pairwise.complete.obs", method = "spearman")
print(round(cor_mat, 3))

cm_long <- as.data.frame(cor_mat) %>%
  tibble::rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "r")

save_plot(
  ggplot(cm_long, aes(var1, var2, fill = r)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
    scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759",
                         limits = c(-1, 1)) +
    labs(x = NULL, y = NULL,
         title = "Corrélation de Spearman entre marques et RNA",
         subtitle = "Une colonne par marque, une ligne par gène") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1)),
  "11_correlation_matrix.pdf", 7, 6)

# ===========================================================================
# 12. EFFET ChIP PAR CLASSE DE RÉGION   (dat_all : c'est le but)
# ===========================================================================
save_plot(
  ggplot(dat_all %>% filter(region_class %in% c("Promoter", "Genic", "Distal")),
         aes(region_class, logFC_chip, fill = region_class)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_boxplot(outlier.size = 0.2) +
    facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_y") +
    scale_fill_brewer(palette = "Pastel1") +
    labs(x = NULL, y = "log2FC ChIP",
         title = "Effet chromatinien par classe de région") +
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "12_chip_by_region_boxplot.pdf")

# ===========================================================================
# RÉCAP
# ===========================================================================
summary_tab <- conc %>%
  dplyr::select(mark, type, n, pct_observed = pct, OR, fisher_p) %>%
  left_join(perm %>% dplyr::select(mark, null_mean, excess, z, p_perm), by = "mark") %>%
  left_join(robust %>% dplyr::select(mark, delta_agg = delta), by = "mark")
if (nrow(auc_tab))
  summary_tab <- summary_tab %>%
  left_join(auc_tab %>% dplyr::select(mark, AUC_signed = AUC), by = "mark")

cat("\n=========================== RÉCAP ===========================\n")
print(as.data.frame(summary_tab %>% mutate(fisher_p = signif(fisher_p, 3))),
      row.names = FALSE)
write_csv(summary_tab, file.path(INT_DIR, "SUMMARY.csv"))

cat("\n==> Figures : ", PLOT_DIR, "\n")
cat("    À regarder EN PREMIER :\n")
cat("      03_concordance_by_abundance.pdf  -> l'effet est-il un artefact ?\n")
cat("      04_concordance_by_region.pdf     -> promoteur ou enhancer ?\n")
cat("      09c_FA_module_scores.pdf         -> NegCtrl doit être plat\n\n")





