#!/usr/bin/env Rscript
# =============================================================================
# 05b_integrate_edgeR.R — Intégration RNA-seq (edgeR, Ensembl) <-> ChIP-seq
# -----------------------------------------------------------------------------
# Adapté à un CSV edgeR dont la 1re colonne (sans nom) est l'identifiant Ensembl.
# Étapes :
#   1. Lit le CSV edgeR (séparateur virgule, 1re colonne = Ensembl).
#   2. Mappe Ensembl -> Entrez + Symbol via org.Hs.eg.db.
#   3. Fusionne avec chaque table ChIP *_gene_level.csv (clé = Entrez/GeneID).
#   4. Produit la table cible + une colonne 'category'.
#   5. Extrait explicitement FXN sur les 4 marques + son expression RNA.
# Sorties : results/integration/integration_<MARK>.csv
#           results/integration/FXN_summary.csv
# =============================================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi); library(dplyr); library(readr)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
ANN_DIR <- file.path(PROJECT_ROOT, "results", "annotation")
INT_DIR <- file.path(PROJECT_ROOT, "results", "integration")
dir.create(INT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Fichier RNA-seq (À VÉRIFIER) ----
RNA_FILE <- "/home/apuhart/Downloads/Bulk-RNAseq/results_test/fibro/control_vs_FA/fibro_control_vs_FA_DESeq2_results.csv"
# Sens attendu par marque (répressif : anticorrélation ; actif : corrélation)
MARK_TYPE <- c(H3K9me3="repressive", H3K27me3="repressive",
               H3K27ac="active", PolII="active")
LFC_TH <- 1; SIG_TH <- 0.05   # seuils indicatifs pour la colonne 'category'

# ---------------------------------------------------------------------------
# 1-2. RNA-seq + conversion Ensembl -> Entrez / Symbol
# ---------------------------------------------------------------------------
rna <- read_csv(RNA_FILE, show_col_types = FALSE)
names(rna)[1] <- "ENSEMBL"                 # 1re colonne sans nom = Ensembl
rna <- rna %>% transmute(
  ENSEMBL = sub("\\..*$", "", ENSEMBL),    # retire un éventuel suffixe .version
  logFC_RNA = as.numeric(log2FoldChange),
  FDR_RNA   = as.numeric(padj),
  PValue_RNA = as.numeric(pvalue))

map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(rna$ENSEMBL),
                             keytype = "ENSEMBL",
                             columns = c("ENTREZID", "SYMBOL")) %>%
  distinct(ENSEMBL, .keep_all = TRUE)
rna <- rna %>% left_join(map, by = "ENSEMBL")
cat(sprintf("RNA : %d gènes, %d mappés vers Entrez.\n",
            nrow(rna), sum(!is.na(rna$ENTREZID))))

# ---------------------------------------------------------------------------
# 3-4. Fusion avec chaque table ChIP gene_level
# ---------------------------------------------------------------------------
classify <- function(d, type) {
  sc <- !is.na(d$FDR_chip) & d$FDR_chip < SIG_TH & abs(d$logFC_chip) >= LFC_TH
  sr <- !is.na(d$FDR_RNA) & d$FDR_RNA < SIG_TH & abs(d$logFC_RNA) >= LFC_TH
  if (type == "repressive") {
    dplyr::case_when(
      sc & sr & d$logFC_chip > 0 & d$logFC_RNA < 0 ~ "gain_mark_down_expr",
      sc & sr & d$logFC_chip < 0 & d$logFC_RNA > 0 ~ "loss_mark_up_expr",
      TRUE ~ "ns_or_other")
  } else {
    dplyr::case_when(
      sc & sr & d$logFC_chip > 0 & d$logFC_RNA > 0 ~ "gain_mark_up_expr",
      sc & sr & d$logFC_chip < 0 & d$logFC_RNA < 0 ~ "loss_mark_down_expr",
      TRUE ~ "ns_or_other")
  }
}

files <- list.files(ANN_DIR, pattern = "_gene_level\\.csv$", full.names = TRUE)
for (f in files) {
  mark <- sub("__.*", "", basename(f))
  chip <- read_csv(f, show_col_types = FALSE) %>%
    transmute(ENTREZID = as.character(GeneID), GeneSymbol,
              signal_chip = mean_signal, logFC_chip = logFC, FDR_chip = FDR)
  
  merged <- rna %>%
    mutate(ENTREZID = as.character(ENTREZID)) %>%
    inner_join(chip, by = "ENTREZID")
  if (nrow(merged) == 0) { cat(sprintf("[%s] aucun gène commun.\n", mark)); next }
  merged$category <- classify(merged, MARK_TYPE[[mark]])
  
  out <- merged %>%
    transmute(GeneSymbol = coalesce(GeneSymbol, SYMBOL), ENSEMBL, ENTREZID,
              logFC_RNA, FDR_RNA,
              !!paste0("signal_", mark) := signal_chip,
              !!paste0("logFC_",  mark) := logFC_chip,
              !!paste0("FDR_",    mark) := FDR_chip,
              category) %>%
    arrange(FDR_RNA)
  write_csv(out, file.path(INT_DIR, sprintf("integration_%s.csv", mark)))
  cat(sprintf("[%s] %d gènes communs -> integration_%s.csv\n",
              mark, nrow(out), mark))
}

# ---------------------------------------------------------------------------
# 5. Focus FXN (toujours, même non significatif)
# ---------------------------------------------------------------------------
FXN_ENS <- "ENSG00000165060"
fxn_rows <- list()
fxn_rna <- rna %>% filter(ENSEMBL == FXN_ENS)
for (f in files) {
  mark <- sub("__.*", "", basename(f))
  chip <- read_csv(f, show_col_types = FALSE)
  row <- chip %>% filter(GeneSymbol == "FXN")
  if (nrow(row)) fxn_rows[[mark]] <- data.frame(
    mark = mark, logFC_chip = row$logFC[1], FDR_chip = row$FDR[1],
    n_regions = row$n_regions[1])
}
fxn_chip <- if (length(fxn_rows)) bind_rows(fxn_rows) else data.frame()
cat("\n===== FXN =====\n")
cat(sprintf("RNA-seq : logFC=%.3f  FDR=%.3f\n",
            ifelse(nrow(fxn_rna), fxn_rna$logFC_RNA[1], NA),
            ifelse(nrow(fxn_rna), fxn_rna$FDR_RNA[1], NA)))
if (nrow(fxn_chip)) { print(fxn_chip, row.names = FALSE)
  write_csv(fxn_chip, file.path(INT_DIR, "FXN_summary.csv"))
} else cat("Aucun pic ChIP assigné à FXN (cf. module ciblé 06/07).\n")

cat("\n==> [05b] Terminé. Tables dans results/integration/\n")