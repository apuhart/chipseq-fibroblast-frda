#!/usr/bin/env Rscript
# =============================================================================
# 16_integration_analysis.R — Intégration RNA-seq x ChIP-seq (FRDA)   [v4]
# -----------------------------------------------------------------------------
# PRÉREQUIS : 04b_recompute_gene_level.R (version corrigée, cf. aggregate_gene)
#
# NOUVEAUTÉS v4
#  [NEW] §0c : contrôle de DÉRIVE GLOBALE. La moyenne génome-entier de
#        logFC_chip est-elle nulle ? Si non, comparer un module à 0 est faux.
#  [NEW] §3  : NULL PAR PERMUTATION *DANS CHAQUE QUINTILE*. Le null global
#        n'est pas 50% ; le null d'un quintile ne l'est pas non plus.
#  [REBUILT] §9c : les modules sont désormais testés en COMPÉTITIF (module vs
#        reste du génome) avec IC95 et correction BH. L'ancienne version
#        comparait une moyenne sur n=1..13 gènes à zéro, sans test ni IC :
#        à n=4, l'IC95 valait +/- 0.49, soit toute la hauteur de l'axe.
#        Les modules à n < MIN_SET_SIZE ne sont plus affichés du tout.
#  [NEW] §9c-bis : sets Reactome via msigdbr (50-200 gènes) — c'est ce qui
#        donne de la puissance. Les listes curées à la main de 5-18 gènes
#        n'en ont aucune une fois intersectées avec ChIP ∩ RNA.
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

CHIP_AGG      <- "mean"   # "mean" | "maxabs"
N_PERM        <- 1000     # permutations, analyse globale
N_PERM_STRAT  <- 500      # permutations, par strate (plus nombreuses -> moins)
FDR_RNA_TH    <- 0.05
LFC_RNA_TH    <- 1
LFC_RNA_MAX   <- 15
MIN_SET_SIZE  <- 8        # en dessous : un "module score" n'a pas de sens
USE_MSIGDB    <- TRUE     # sets Reactome (recommandé : 50-200 gènes)

PAL <- c(repressive = "#4E79A7", active = "#E15759")

save_plot <- function(p, filename, width = 11, height = 8) {
  ggsave(file.path(PLOT_DIR, filename), p, width = width, height = height)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Helpers de concordance — TYPE-AWARE
# ---------------------------------------------------------------------------
# répressive : concordant = signes OPPOSÉS | active : concordant = MÊME signe
concord_class <- function(lfc_chip, lfc_rna, type) {
  p <- lfc_chip * lfc_rna
  ifelse(is.na(p) | p == 0, NA_character_,
         ifelse((type == "active" & p > 0) | (type == "repressive" & p < 0),
                "Concordant", "Discordant"))
}
pct_concordant <- function(lfc_chip, lfc_rna, type) {
  cc <- concord_class(lfc_chip, lfc_rna, type); cc <- cc[!is.na(cc)]
  if (!length(cc)) return(NA_real_)
  100 * mean(cc == "Concordant")
}
# Null par permutation : mélange les étiquettes, préserve les deux marges.
perm_null <- function(lfc_chip, lfc_rna, type, n_perm) {
  obs <- pct_concordant(lfc_chip, lfc_rna, type)
  nv  <- replicate(n_perm, pct_concordant(sample(lfc_chip), lfc_rna, type))
  list(observed = obs, null_mean = mean(nv), null_sd = sd(nv),
       excess = obs - mean(nv), z = (obs - mean(nv)) / sd(nv),
       p = (1 + sum(abs(nv - mean(nv)) >= abs(obs - mean(nv)))) / (n_perm + 1))
}

# ---------------------------------------------------------------------------
# 0a. RNA-seq
# ---------------------------------------------------------------------------
cat("\n=========================== CHARGEMENT ===========================\n")
cat(sprintf("Agrégation ChIP : %s | MIN_SET_SIZE : %d\n", CHIP_AGG, MIN_SET_SIZE))

rna_raw <- read_csv(RNA_FILE, show_col_types = FALSE)
names(rna_raw)[1] <- "GENE_ID"
need <- c("log2FoldChange", "pvalue", "padj")
if (length(setdiff(need, names(rna_raw))))
  stop("Colonnes absentes du CSV RNA : ", paste(setdiff(need, names(rna_raw)), collapse = ", "))

rna <- rna_raw %>%
  transmute(GENE_ID    = sub("\\..*$", "", as.character(GENE_ID)),
            baseMean   = if ("baseMean" %in% names(rna_raw)) as.numeric(baseMean) else NA_real_,
            logFC_RNA  = as.numeric(log2FoldChange),
            FDR_RNA    = as.numeric(padj),
            PValue_RNA = as.numeric(pvalue)) %>%
  filter(!is.na(GENE_ID), GENE_ID != "")

n_extreme <- sum(abs(rna$logFC_RNA) > LFC_RNA_MAX, na.rm = TRUE)
if (n_extreme) {
  cat(sprintf("!! %d gènes |logFC_RNA| > %d écartés (comptage nul dans un groupe).\n",
              n_extreme, LFC_RNA_MAX))
  rna <- rna %>% filter(is.na(logFC_RNA) | abs(logFC_RNA) <= LFC_RNA_MAX)
}

is_ens <- grepl("^ENSG[0-9]+$", rna$GENE_ID)
safe_select <- function(k, kt, cols) {
  ok <- intersect(unique(k), AnnotationDbi::keys(org.Hs.eg.db, kt))
  if (!length(ok)) return(NULL)
  suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = ok, keytype = kt, columns = cols))
}
m_sym <- safe_select(rna$GENE_ID[!is_ens], "SYMBOL", c("ENTREZID", "ENSEMBL"))
m_sym <- if (is.null(m_sym)) tibble() else m_sym %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% transmute(GENE_ID = SYMBOL, ENTREZID, SYMBOL)
m_ens <- safe_select(rna$GENE_ID[is_ens], "ENSEMBL", c("ENTREZID", "SYMBOL"))
m_ens <- if (is.null(m_ens)) tibble() else m_ens %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>% transmute(GENE_ID = ENSEMBL, ENTREZID, SYMBOL)

rna <- rna %>%
  left_join(bind_rows(m_sym, m_ens) %>% distinct(GENE_ID, .keep_all = TRUE), by = "GENE_ID") %>%
  filter(!is.na(ENTREZID)) %>% mutate(ENTREZID = as.character(ENTREZID))

n_before <- nrow(rna)
rna <- rna %>% arrange(FDR_RNA, PValue_RNA) %>% distinct(ENTREZID, .keep_all = TRUE)
cat(sprintf("RNA  : %d gènes (%d doublons Entrez écartés) | %d à FDR<%.2f\n",
            nrow(rna), n_before - nrow(rna),
            sum(rna$FDR_RNA < FDR_RNA_TH, na.rm = TRUE), FDR_RNA_TH))

# ---------------------------------------------------------------------------
# 0b. ChIP + jointure
# ---------------------------------------------------------------------------
files <- list.files(ANN_DIR, pattern = "_gene_level_fixed\\.csv$", full.names = TRUE)
if (!length(files)) stop("Aucun *_gene_level_fixed.csv -> lance 04b d'abord.")

chip_raw <- lapply(files, function(f) {
  tag <- sub("_gene_level_fixed$", "", tools::file_path_sans_ext(basename(f)))
  read_csv(f, show_col_types = FALSE) %>%
    mutate(mark = sub("__.*$", "", tag), ENTREZID = as.character(GeneID))
}) %>% bind_rows() %>% distinct(mark, ENTREZID, region_class, .keep_all = TRUE)

if (!"mean_logFC" %in% names(chip_raw))
  stop("Colonne mean_logFC absente -> tes tables viennent du 04b NON corrigé.\n",
       "Relance 04b_recompute_gene_level.R (version avec aggregate_gene corrigée).")

# Sanity check : mean_logFC doit DIFFÉRER de logFC pour les gènes multi-régions.
# Si elles sont identiques partout, le 04b bugué (masquage dplyr) a été utilisé.
.multi <- chip_raw %>% filter(region_class == "All", n_regions > 1)
if (nrow(.multi) > 0) {
  .id <- mean(abs(.multi$mean_logFC - .multi$logFC) < 1e-12, na.rm = TRUE)
  if (.id > 0.99)
    stop("mean_logFC == logFC pour ", round(100 * .id, 1), "% des gènes multi-régions.\n",
         "-> tes tables viennent du 04b AVANT correction du masquage dplyr.\n",
         "   Relance 04b_recompute_gene_level.R.")
  cat(sprintf("ChIP : %d gènes multi-régions, mean_logFC != logFC pour %.1f%% -> 04b OK.\n",
              nrow(.multi), 100 * (1 - .id)))
}

build_dat <- function(chip_df, agg) {
  chip_df %>%
    mutate(logFC_chip = if (agg == "mean") mean_logFC else logFC) %>%
    inner_join(rna, by = "ENTREZID", relationship = "many-to-one") %>%
    transmute(mark, type = unname(MARK_TYPE[mark]), region_class,
              gene = coalesce(GeneSymbol, SYMBOL), ENTREZID, n_regions, mean_signal,
              logFC_chip, FDR_chip = FDR, baseMean, logFC_RNA, PValue_RNA, FDR_RNA) %>%
    filter(!is.na(logFC_chip), !is.na(logFC_RNA))
}

dat_all <- build_dat(chip_raw, CHIP_AGG)
dat     <- dat_all %>% filter(region_class == "All")   # 1 ligne / gène / marque

# ---------------------------------------------------------------------------
# 0c. DIAGNOSTIC + DÉRIVE GLOBALE
# ---------------------------------------------------------------------------
cat("\n=========================== DIAGNOSTIC ===========================\n")
diag <- dat %>%
  group_by(mark, type) %>%
  summarise(n_genes = n(),
            max_abs_lfc_chip = round(max(abs(logFC_chip), na.rm = TRUE), 3),
            sig_chip_05 = sum(FDR_chip < 0.05, na.rm = TRUE),
            sig_rna_05  = sum(FDR_RNA < 0.05, na.rm = TRUE),
            .groups = "drop")
print(as.data.frame(diag), row.names = FALSE)

if (any(diag$max_abs_lfc_chip < 0.1))
  stop("\n!! logFC ChIP ÉCRASÉ : relance 04b_recompute_gene_level.R\n")

# DÉRIVE GLOBALE : si la moyenne génome-entier de logFC_chip n'est pas ~0,
# tout module comparé à ZÉRO montrera un effet qui n'est que cette dérive.
# C'est pourquoi §9c teste en COMPÉTITIF (module vs reste), pas contre 0.
cat("\n--- Dérive globale du logFC ChIP ---\n")
drift <- dat %>% group_by(mark, type) %>%
  summarise(mean_all = round(mean(logFC_chip, na.rm = TRUE), 4),
            med_all  = round(median(logFC_chip, na.rm = TRUE), 4),
            sd_all   = round(sd(logFC_chip, na.rm = TRUE), 3),
            pct_neg  = round(100 * mean(logFC_chip < 0, na.rm = TRUE), 1),
            .groups = "drop")
print(as.data.frame(drift), row.names = FALSE)
cat("mean_all loin de 0 = dérive globale -> ne JAMAIS comparer un module à 0.\n")
cat("pct_neg loin de 50 = marges asymétriques -> le null de concordance != 50%.\n")
write_csv(drift, file.path(INT_DIR, "global_drift.csv"))

# ---------------------------------------------------------------------------
# 0d. ROBUSTESSE DE L'AGRÉGATION
# ---------------------------------------------------------------------------
cat("\n=================== 0d. ROBUSTESSE DE L'AGRÉGATION ===================\n")
robust <- lapply(c("mean", "maxabs"), function(agg) {
  build_dat(chip_raw, agg) %>% filter(region_class == "All") %>%
    group_by(mark, type) %>%
    summarise(pct = round(pct_concordant(logFC_chip, logFC_RNA, type[1]), 1),
              .groups = "drop") %>% mutate(agg = agg)
}) %>% bind_rows() %>%
  pivot_wider(names_from = agg, values_from = pct) %>%
  mutate(delta = round(mean - maxabs, 1))
print(as.data.frame(robust), row.names = FALSE)
cat("\n|delta| < 2 -> robuste. |delta| > 5 -> dépend de la région choisie.\n")
if (all(abs(robust$delta) < 1e-9))
  warning("delta == 0 pour TOUTES les marques : suspect. Vérifie que 04b a bien\n",
          "été relancé après correction (mean_logFC ne doit pas copier logFC).")
write_csv(robust, file.path(INT_DIR, "aggregation_robustness.csv"))

# ===========================================================================
# 1. CONCORDANCE GLOBALE
# ===========================================================================
# Le % brut dépend des marges ; l'ODDS RATIO n'en dépend pas. Quand les deux
# divergent (cas de PolII), croire l'OR.
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
cat("\nOR attendu : >1 si ACTIVE, <1 si RÉPRESSIVE.\n")

# ===========================================================================
# 2. NULL PAR PERMUTATION (global)
# ===========================================================================
cat(sprintf("\n=================== 2. NULL PAR PERMUTATION (%d) ===================\n", N_PERM))
perm <- lapply(split(dat, dat$mark), function(d) {
  d <- d %>% filter(logFC_chip != 0, logFC_RNA != 0)
  r <- perm_null(d$logFC_chip, d$logFC_RNA, d$type[1], N_PERM)
  tibble(mark = d$mark[1], type = d$type[1],
         observed = round(r$observed, 2), null_mean = round(r$null_mean, 2),
         excess = round(r$excess, 2), null_sd = round(r$null_sd, 3),
         z = round(r$z, 1), p_perm = r$p)
}) %>% bind_rows() %>% arrange(desc(excess))
print(as.data.frame(perm), row.names = FALSE)
write_csv(perm, file.path(INT_DIR, "permutation_null.csv"))
cat("\nexcess = TAILLE de l'effet (points) | z = SIGNIFICATIVITÉ.\n")

perm_plot <- lapply(split(dat, dat$mark), function(d) {
  d <- d %>% filter(logFC_chip != 0, logFC_RNA != 0); ty <- d$type[1]
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
         subtitle = "Le null n'est PAS à 50% : il dépend des marges de chaque marque") +
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "02_permutation_null.pdf")

# ===========================================================================
# 3. STRATIFICATION PAR ABONDANCE + NULL PAR QUINTILE  <<< CONTRÔLE DÉCISIF
# ===========================================================================
# Confondant redouté : les gènes fortement couverts ont plus de marque active
# ET plus de puissance RNA -> la concordance pourrait n'être que ça.
# SIGNATURE D'UN ARTEFACT : Q1 ~ null, Q5 >> null (pente MONTANTE).
# SIGNATURE D'UN VRAI EFFET : excès homogène sur les 5 quintiles.
#
# Le null étant marge-dépendant, il est recalculé DANS CHAQUE QUINTILE.
cat(sprintf("\n=========== 3. STRATIFICATION PAR ABONDANCE (+ null x %d) ===========\n",
            N_PERM_STRAT))

strat <- lapply(split(dat, dat$mark), function(d) {
  ty <- d$type[1]
  d <- d %>% mutate(q_signal = ntile(mean_signal, 5))
  lapply(1:5, function(q) {
    dq <- d %>% filter(q_signal == q, logFC_chip != 0, logFC_RNA != 0)
    if (nrow(dq) < 30) return(NULL)
    r <- perm_null(dq$logFC_chip, dq$logFC_RNA, ty, N_PERM_STRAT)
    tibble(mark = d$mark[1], type = ty, q_signal = q, n = nrow(dq),
           observed = round(r$observed, 1), null_mean = round(r$null_mean, 1),
           excess = round(r$excess, 1), z = round(r$z, 1), p_perm = r$p)
  }) %>% bind_rows()
}) %>% bind_rows()

print(as.data.frame(strat), row.names = FALSE)
write_csv(strat, file.path(INT_DIR, "concordance_by_abundance.csv"))

# Test d'homogénéité de l'excès entre quintiles
cat("\n--- Homogénéité de l'excès entre quintiles ---\n")
hom <- strat %>% group_by(mark, type) %>%
  summarise(excess_min = min(excess), excess_max = max(excess),
            excess_range = round(max(excess) - min(excess), 1),
            slope_Q1_Q5 = round(excess[q_signal == 5] - excess[q_signal == 1], 1),
            .groups = "drop") %>%
  mutate(verdict = case_when(
    slope_Q1_Q5 > 4  ~ "PENTE MONTANTE -> suspect (artefact d'abondance ?)",
    excess_range < 4 ~ "homogène -> effet réel",
    TRUE             ~ "hétérogène mais pente non montante -> acceptable"))
print(as.data.frame(hom), row.names = FALSE)

save_plot(
  ggplot(strat, aes(factor(q_signal))) +
    geom_col(aes(y = observed, fill = type), width = 0.7, alpha = 0.9) +
    geom_point(aes(y = null_mean), shape = 95, size = 8, colour = "black") +
    geom_text(aes(y = observed, label = sprintf("%.0f%%\nz=%.0f", observed, z)),
              vjust = -0.25, size = 2.4) +
    scale_fill_manual(values = PAL) +
    facet_wrap(~ paste0(mark, " (", type, ")")) +
    coord_cartesian(ylim = c(0, 80)) +
    labs(x = "Quintile de signal ChIP (1 = faible couverture → 5 = forte)",
         y = "% concordant",
         title = "La concordance tient-elle à toutes les abondances ?",
         subtitle = "Trait noir = null par permutation DU QUINTILE (pas 50%)") +
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "03_concordance_by_abundance.pdf", 10, 7)

# ---- 3b. Concordance par décile RNA ---------------------------------------
binned_conc <- dat %>%
  group_by(mark) %>% mutate(bin = ntile(logFC_RNA, 10)) %>%
  group_by(mark, type, bin) %>%
  summarise(n = n(), pct = round(pct_concordant(logFC_chip, logFC_RNA, type[1]), 1),
            mean_rna = round(mean(logFC_RNA, na.rm = TRUE), 3),
            mean_chip = round(mean(logFC_chip, na.rm = TRUE), 3), .groups = "drop")
write_csv(binned_conc, file.path(INT_DIR, "concordance_by_rna_decile.csv"))

save_plot(
  ggplot(binned_conc, aes(bin, pct, colour = type)) +
    geom_hline(yintercept = 50, linetype = 2, colour = "grey60") +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    scale_x_continuous(breaks = 1:10) + scale_colour_manual(values = PAL) +
    facet_wrap(~ paste0(mark, " (", type, ")")) +
    labs(x = "Décile RNA (1 = plus réprimé → 10 = plus induit)", y = "% concordant",
         title = "Concordance par décile RNA",
         subtitle = "Se renforce-t-elle aux extrêmes ?") +
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "03b_concordance_by_rna_decile.pdf")

# ===========================================================================
# 4. PROMOTEUR vs DISTAL   (dat_all : c'est le but de la section)
# ===========================================================================
cat("\n=================== 4. PROMOTEUR vs DISTAL ===================\n")
by_reg <- lapply(split(dat_all %>% filter(region_class %in% c("Promoter","Genic","Distal")),
                       ~ mark + region_class), function(d) {
                         if (nrow(d) < 50) return(NULL)
                         d <- d %>% filter(logFC_chip != 0, logFC_RNA != 0)
                         if (nrow(d) < 50) return(NULL)
                         r <- perm_null(d$logFC_chip, d$logFC_RNA, d$type[1], N_PERM_STRAT)
                         tibble(mark = d$mark[1], type = d$type[1], region_class = d$region_class[1],
                                n = nrow(d), observed = round(r$observed, 1),
                                null_mean = round(r$null_mean, 1), excess = round(r$excess, 1),
                                z = round(r$z, 1), p_perm = r$p)
                       }) %>% bind_rows()

if (nrow(by_reg)) {
  by_reg <- by_reg %>% arrange(mark, region_class)
  print(as.data.frame(by_reg), row.names = FALSE)
  write_csv(by_reg, file.path(INT_DIR, "concordance_by_region_class.csv"))
  save_plot(
    ggplot(by_reg, aes(region_class)) +
      geom_col(aes(y = observed, fill = type), width = 0.7) +
      geom_point(aes(y = null_mean), shape = 95, size = 8, colour = "black") +
      geom_text(aes(y = observed, label = sprintf("%.0f%%\nexcès %+.1f\nn=%d",
                                                  observed, excess, n)),
                vjust = -0.2, size = 2.3) +
      scale_fill_manual(values = PAL) +
      facet_wrap(~ paste0(mark, " (", type, ")")) +
      coord_cartesian(ylim = c(0, 90)) +
      labs(x = NULL, y = "% concordant",
           title = "Concordance par classe de région",
           subtitle = "Trait noir = null du groupe | Promoteur = direct, Distal = enhancer") +
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
            se_chip = sd(logFC_chip, na.rm = TRUE) / sqrt(n()), n = n(), .groups = "drop") %>%
  left_join(drift %>% dplyr::select(mark, mean_all), by = "mark")

save_plot(
  ggplot(binned, aes(bin, mean_chip)) +
    geom_hline(aes(yintercept = mean_all), linetype = 2, colour = "grey40") +
    geom_errorbar(aes(ymin = mean_chip - 1.96 * se_chip,
                      ymax = mean_chip + 1.96 * se_chip), width = 0.2, colour = "grey50") +
    geom_line(aes(colour = type), linewidth = 0.8) +
    geom_point(aes(colour = type), size = 2.5) +
    scale_x_continuous(breaks = 1:10) + scale_colour_manual(values = PAL) +
    facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_y") +
    labs(x = "Décile de log2FC RNA (1 = plus réprimé → 10 = plus induit)",
         y = "log2FC ChIP moyen (± IC95)",
         title = "Le signal chromatinien suit-il l'expression ?",
         subtitle = "Trait pointillé = moyenne génome-entier (pas zéro)") +
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "05_binned_trend.pdf")

# ===========================================================================
# 6. SCATTER
# ===========================================================================
lab <- dat %>% group_by(mark) %>%
  filter(gene == "FXN" | rank(FDR_RNA, ties.method = "first") <= 8) %>% ungroup()

save_plot(
  ggplot(dat, aes(logFC_chip, logFC_RNA)) +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(aes(colour = FDR_RNA < FDR_RNA_TH), alpha = 0.35, size = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#E15759", linewidth = 0.7) +
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
# Question : parmi les gènes RNA-significatifs, le logFC_chip SIGNÉ prédit-il
# le SENS (up/down) ? (Et non : |logFC_chip| prédit-il "significatif" — ce qui
# mesurerait surtout la qualité de mesure.)
auc_tab <- tibble()
if (requireNamespace("pROC", quietly = TRUE)) {
  suppressPackageStartupMessages(library(pROC))
  cat("\n=================== 8. ROC SIGNÉE (direction) ===================\n")
  roc_dir <- file.path(PLOT_DIR, "ROC"); dir.create(roc_dir, showWarnings = FALSE)
  auc_tab <- lapply(split(dat, dat$mark), function(d) {
    dd <- d %>% filter(FDR_RNA < FDR_RNA_TH, abs(logFC_RNA) > LFC_RNA_TH)
    if (nrow(dd) < 20 || n_distinct(sign(dd$logFC_RNA)) < 2) {
      cat(sprintf("  [%s] n=%d -> trop peu pour une ROC.\n", d$mark[1], nrow(dd)))
      return(NULL)
    }
    resp <- factor(ifelse(dd$logFC_RNA > 0, "Up", "Down"), levels = c("Down", "Up"))
    pred <- if (dd$type[1] == "repressive") -dd$logFC_chip else dd$logFC_chip
    r <- pROC::roc(resp, pred, quiet = TRUE, direction = "<")
    pdf(file.path(roc_dir, paste0("ROC_", dd$mark[1], ".pdf")), width = 5, height = 5)
    plot(r, legacy.axes = TRUE, col = PAL[[dd$type[1]]], lwd = 3,
         main = sprintf("%s (%s)\nAUC = %.3f  |  n = %d", dd$mark[1], dd$type[1],
                        as.numeric(pROC::auc(r)), nrow(dd)))
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
        scale_fill_manual(values = PAL) + coord_flip(ylim = c(0.3, 1)) +
        labs(x = NULL, y = "ROC AUC (± IC95)",
             title = "La chromatine prédit-elle le SENS du changement RNA ?") +
        theme_bw(base_size = 11),
      "08_AUC_signed.pdf", 7, 4.5)
  }
} else cat("\n[skip] pROC absent.\n")

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
      ego <- tryCatch(enrichGO(gene = ids, universe = universe_e, OrgDb = org.Hs.eg.db,
                               keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                               pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE),
                      error = function(e) NULL)
      n_term <- if (!is.null(ego)) nrow(as.data.frame(ego)) else 0
      cat(sprintf("  -> %d termes\n", n_term))
      if (n_term > 0) go_res[[g]] <- as.data.frame(ego) %>% mutate(group = g)
    }
  } else cat("[skip] clusterProfiler absent.\n")
}

# GeneRatio est une chaîne "5/100" -> à parser, sinon ggplot en fait un facteur.
if (length(go_res)) {
  go_df <- bind_rows(go_res) %>%
    mutate(GeneRatio_num = vapply(strsplit(GeneRatio, "/"),
                                  function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1)))
  write_csv(go_df, file.path(INT_DIR, sprintf("%s_GO_by_group.csv", FOCUS_MARK)))
  for (g in unique(go_df$group)) {
    sub <- go_df %>% filter(group == g) %>% arrange(p.adjust) %>% slice_head(n = 20)
    save_plot(
      ggplot(sub, aes(GeneRatio_num, reorder(Description, GeneRatio_num),
                      size = Count, colour = p.adjust)) +
        geom_point() + scale_colour_gradient(low = "#E15759", high = "#4E79A7") +
        labs(x = "Gene ratio", y = NULL, title = paste("GO:BP —", FOCUS_MARK, "—", g),
             subtitle = "Top 20 termes") +
        theme_bw(base_size = 11) + theme(axis.text.y = element_text(size = 7)),
      sprintf("09b_GO_%s_%s.pdf", FOCUS_MARK, g), 9, 6)
  }
}

# ===========================================================================
# 9c. GENE SETS — TEST COMPÉTITIF   [SECTION ENTIÈREMENT REFAITE]
# ===========================================================================
# POURQUOI L'ANCIENNE VERSION ÉTAIT INUTILISABLE :
#   - moyennes calculées sur n = 1 à 13 gènes, sans IC ni test
#     (à n=4 l'IC95 vaut +/- 0.49, soit toute la hauteur de l'axe)
#   - comparaison à ZÉRO alors que la moyenne génome-entier ne l'est pas
#     forcément (cf. §0c) -> on mesurait la dérive globale, pas le module
#
# CE QUI EST FAIT ICI :
#   - test COMPÉTITIF : module vs RESTE du génome (Wilcoxon rang-somme)
#   - IC95 sur la moyenne du set, référence = moyenne du reste
#   - correction BH sur l'ensemble des tests
#   - les sets à n < MIN_SET_SIZE ne sont PAS affichés
#   - sets Reactome (50-200 gènes) en plus des listes curées : c'est ce qui
#     donne la puissance. Une liste de 5 gènes n'en aura jamais.
cat("\n=================== 9c. GENE SETS (test compétitif) ===================\n")

curated_sets <- list(
  FeS_assembly = c("FXN","NFS1","ISCU","ISCA1","ISCA2","FDX1","FDX2","NUBPL","NFU1",
                   "BOLA1","BOLA3","LYRM4","GLRX5","LIAS","LIPT1","LIPT2","HSPA9","HSCB"),
  Iron_homeostasis = c("TFRC","SLC40A1","FTH1","FTL","SLC11A2","ACO1","ACO2","IREB2"),
  OXPHOS = c("NDUFS1","NDUFS2","NDUFA9","NDUFB8","SDHA","SDHB","UQCRC1","UQCRC2",
             "COX4I1","COX5A","ATP5F1A"),
  Oxidative_stress = c("NFE2L2","HMOX1","SOD1","SOD2","GPX1","GPX4","NQO1","GCLC",
                       "CAT","PRDX1","PRDX3","TXN","TXNRD1"),
  Mito_biogenesis = c("PPARGC1A","TFAM","NRF1","ESRRA","PPRC1")
)

gene_sets <- curated_sets
if (USE_MSIGDB && requireNamespace("msigdbr", quietly = TRUE)) {
  msig <- tryCatch({
    m <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
    keep <- grepl("IRON|SULFUR|RESPIRATORY_ELECTRON|TCA_CYCLE|MITOCHONDRIAL|OXIDATIVE_STRESS|DETOXIFICATION",
                  m$gs_name)
    split(unique(m$gene_symbol[keep]), m$gs_name[keep])
  }, error = function(e) NULL)
  if (!is.null(msig) && length(msig)) {
    msig <- msig[lengths(msig) >= 20 & lengths(msig) <= 300]
    gene_sets <- c(gene_sets, msig)
    cat(sprintf("  + %d sets Reactome (20-300 gènes) chargés.\n", length(msig)))
  } else cat("  [msigdbr indisponible ou aucun set retenu]\n")
}

set_test <- lapply(split(dat, dat$mark), function(d) {
  lapply(names(gene_sets), function(sn) {
    inset <- d$gene %in% gene_sets[[sn]]
    n_in <- sum(inset)
    if (n_in < MIN_SET_SIZE) return(NULL)
    w <- suppressWarnings(wilcox.test(d$logFC_chip[inset], d$logFC_chip[!inset]))
    tibble(mark = d$mark[1], type = d$type[1], set = sn, n = n_in,
           mean_set  = mean(d$logFC_chip[inset], na.rm = TRUE),
           se_set    = sd(d$logFC_chip[inset], na.rm = TRUE) / sqrt(n_in),
           mean_rest = mean(d$logFC_chip[!inset], na.rm = TRUE),   # LE VRAI NULL
           delta     = mean(d$logFC_chip[inset], na.rm = TRUE) -
             mean(d$logFC_chip[!inset], na.rm = TRUE),
           pct_conc  = round(pct_concordant(d$logFC_chip[inset],
                                            d$logFC_RNA[inset], d$type[1]), 1),
           p_wilcox  = w$p.value)
  }) %>% bind_rows()
}) %>% bind_rows()

if (nrow(set_test)) {
  set_test <- set_test %>%
    mutate(p_adj = p.adjust(p_wilcox, "BH")) %>%
    arrange(p_adj) %>%
    mutate(across(c(mean_set, se_set, mean_rest, delta), ~ round(.x, 4)))
  print(as.data.frame(set_test %>%
                        mutate(p_wilcox = signif(p_wilcox, 3), p_adj = signif(p_adj, 3)) %>%
                        dplyr::select(mark, set, n, mean_set, mean_rest, delta, pct_conc, p_adj) %>%
                        slice_head(n = 25)), row.names = FALSE)
  write_csv(set_test, file.path(INT_DIR, "geneset_competitive_test.csv"))
  
  cat(sprintf("\n%d sets testés (n >= %d) | %d à BH < 0.05\n",
              nrow(set_test), MIN_SET_SIZE, sum(set_test$p_adj < 0.05)))
  cat("delta = moyenne(set) - moyenne(reste). C'est CE chiffre qui compte,\n")
  cat("pas la moyenne du set comparée à zéro.\n")
  
  cur_only <- set_test %>% filter(set %in% names(curated_sets))
  if (nrow(cur_only)) {
    save_plot(
      ggplot(cur_only, aes(set, mean_set, fill = type)) +
        geom_col(width = 0.65) +
        geom_errorbar(aes(ymin = mean_set - 1.96 * se_set,
                          ymax = mean_set + 1.96 * se_set), width = 0.18) +
        geom_hline(aes(yintercept = mean_rest), linetype = 2, colour = "black") +
        geom_text(aes(label = sprintf("n=%d\nBH=%.2g", n, p_adj)),
                  vjust = -0.4, size = 2.2) +
        scale_fill_manual(values = PAL) +
        facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_y") +
        labs(x = NULL, y = "log2FC ChIP moyen (± IC95)",
             title = "Modules FRDA — test compétitif vs reste du génome",
             subtitle = "Trait noir = moyenne du RESTE des gènes (le vrai null, pas zéro)") +
        theme_bw(base_size = 11) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 40, hjust = 1, size = 7)),
      "09c_geneset_competitive.pdf")
  }
  
  top_sets <- set_test %>% filter(p_adj < 0.05) %>% slice_head(n = 20)
  if (nrow(top_sets))
    save_plot(
      ggplot(top_sets, aes(delta, reorder(str_trunc(set, 45), delta), fill = type)) +
        geom_col() + geom_vline(xintercept = 0, linetype = 2) +
        scale_fill_manual(values = PAL) +
        labs(x = "delta = moyenne(set) - moyenne(reste)", y = NULL,
             title = "Sets significativement décalés (BH < 0.05)") +
        theme_bw(base_size = 11) + theme(axis.text.y = element_text(size = 7)),
      "09c2_geneset_significant.pdf", 9, 6)
} else cat("Aucun set n'atteint MIN_SET_SIZE = ", MIN_SET_SIZE, ".\n", sep = "")

# ===========================================================================
# 10. FXN
# ===========================================================================
cat("\n=========================== 10. FXN ===========================\n")
fxn <- dat_all %>% filter(gene == "FXN") %>%
  dplyr::select(mark, type, region_class, n_regions, logFC_chip, FDR_chip,
                logFC_RNA, FDR_RNA)
if (nrow(fxn)) {
  print(as.data.frame(fxn), row.names = FALSE)
  write_csv(fxn, file.path(INT_DIR, "FXN_cross.csv"))
  fxn_one <- dat %>% filter(gene == "FXN")
  if (nrow(fxn_one))
    save_plot(
      ggplot(fxn_one, aes(mark, logFC_chip, fill = type)) +
        geom_col(width = 0.65) +
        geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
        geom_text(aes(label = sprintf("%.3f", logFC_chip)), vjust = -0.3, size = 3) +
        scale_fill_manual(values = PAL) +
        labs(x = NULL, y = "log2FC ChIP",
             title = "Locus FXN à travers les marques",
             subtitle = sprintf("RNA log2FC = %.2f (FDR = %.3g) | n_regions = %s",
                                fxn_one$logFC_RNA[1], fxn_one$FDR_RNA[1],
                                paste(fxn_one$n_regions, collapse = "/"))) +
        theme_bw(base_size = 11) + theme(legend.position = "none"),
      "10b_FXN_summary.pdf", 8, 4.5)
  cat("\nFXN n'a de pic que sur les marques narrow. Pour les broad au locus FXN\n",
      "-> module ciblé 06/07 (l'assignation pic->gène ne convient pas aux domaines).\n")
} else cat("FXN absent des tables -> module ciblé 06/07.\n")

fes <- dat %>% filter(gene %in% curated_sets$FeS_assembly) %>%
  dplyr::select(gene, mark, type, logFC_RNA, FDR_RNA, logFC_chip)
if (nrow(fes)) {
  write_csv(fes, file.path(INT_DIR, "FeS_module_cross.csv"))
  save_plot(
    ggplot(fes, aes(reorder(gene, logFC_chip), logFC_chip, fill = type)) +
      geom_col(width = 0.65) + geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
      scale_fill_manual(values = PAL) +
      facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_x") +
      labs(x = NULL, y = "log2FC ChIP", title = "Module d'assemblage Fe-S",
           subtitle = "Cœur de la cascade frataxine") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "none"),
    "10c_FeS_module.pdf", 12, 7)
}

# ===========================================================================
# 11. CORRÉLATION ENTRE MARQUES
# ===========================================================================
cat("\n=================== 11. CORRÉLATION ENTRE MARQUES ===================\n")
wide <- dat %>% dplyr::select(ENTREZID, mark, logFC_chip) %>%
  pivot_wider(names_from = mark, values_from = logFC_chip) %>%
  left_join(dat %>% distinct(ENTREZID, logFC_RNA), by = "ENTREZID")
cor_mat <- wide %>% dplyr::select(-ENTREZID) %>%
  cor(use = "pairwise.complete.obs", method = "spearman")
print(round(cor_mat, 3))

cm_long <- as.data.frame(cor_mat) %>% tibble::rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "r")
save_plot(
  ggplot(cm_long, aes(var1, var2, fill = r)) +
    geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
    scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759", limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, title = "Corrélation de Spearman (marques x RNA)",
         subtitle = "Une colonne par marque, une ligne par gène") +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 40, hjust = 1)),
  "11_correlation_matrix.pdf", 7, 6)

# ===========================================================================
# 12. EFFET ChIP PAR CLASSE DE RÉGION
# ===========================================================================
save_plot(
  ggplot(dat_all %>% filter(region_class %in% c("Promoter", "Genic", "Distal")),
         aes(region_class, logFC_chip, fill = region_class)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_boxplot(outlier.size = 0.2) +
    facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_y") +
    scale_fill_brewer(palette = "Pastel1") +
    labs(x = NULL, y = "log2FC ChIP", title = "Effet chromatinien par classe de région") +
    theme_bw(base_size = 11) + theme(legend.position = "none"),
  "12_chip_by_region_boxplot.pdf")

# ===========================================================================
# RÉCAP
# ===========================================================================
summary_tab <- conc %>%
  dplyr::select(mark, type, n, pct_observed = pct, OR, fisher_p) %>%
  left_join(perm %>% dplyr::select(mark, null_mean, excess, z, p_perm), by = "mark") %>%
  left_join(robust %>% dplyr::select(mark, delta_agg = delta), by = "mark") %>%
  left_join(hom %>% dplyr::select(mark, strat_verdict = verdict), by = "mark") %>%
  left_join(drift %>% dplyr::select(mark, global_drift = mean_all), by = "mark")
if (nrow(auc_tab))
  summary_tab <- summary_tab %>% left_join(auc_tab %>% dplyr::select(mark, AUC = AUC), by = "mark")

cat("\n=========================== RÉCAP ===========================\n")
print(as.data.frame(summary_tab %>% mutate(fisher_p = signif(fisher_p, 3))), row.names = FALSE)
write_csv(summary_tab, file.path(INT_DIR, "SUMMARY.csv"))

cat("\n==> Figures : ", PLOT_DIR, "\n")
cat("    ORDRE DE LECTURE :\n")
cat("      03_concordance_by_abundance.pdf   -> l'effet est-il un artefact ?\n")
cat("      04_concordance_by_region.pdf      -> promoteur ou enhancer ?\n")
cat("      02_permutation_null.pdf           -> le contrôle qui valide le reste\n")
cat("      09c_geneset_competitive.pdf       -> modules, avec IC et vrai null\n")
cat("      05_binned_trend.pdf / 06_scatter.pdf  -> illustration\n\n")
