#!/usr/bin/env Rscript
# =============================================================================
# 05_integrate_rnaseq.R — Fusion ChIP-seq <-> RNA-seq
# -----------------------------------------------------------------------------
# Fusionne, pour chaque marque/contraste, la table niveau gène ChIP avec tes
# DEGs RNA-seq, et classe les gènes selon le sens attendu :
#   - marque RÉPRESSIVE (H3K9me3, H3K27me3) : anticorrélation attendue
#       (gain de marque + baisse d'expression ; perte de marque + hausse).
#   - marque ACTIVE (H3K27ac, PolII) : corrélation positive attendue.
#
# Produit la table cible :
#   GeneSymbol | logFC_RNA | padj_RNA | logFC_<MARK> | FDR_<MARK> | category
#
# (Optionnel) Détection de "restauration" si deux contrastes existent pour la
# même marque (ex: FRDA-vs-Control ET Edited-vs-Unedited) : voir bloc en bas.
#
# Lancer :  Rscript scripts/05_integrate_rnaseq.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(purrr)
})

# ---------------------------------------------------------------------------
# CONFIG — À ADAPTER À TON FICHIER RNA-seq
# ---------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
RES     <- file.path(PROJECT_ROOT, "results")
ANN_DIR <- file.path(RES, "annotation")
INT_DIR <- file.path(RES, "integration")
dir.create(INT_DIR, showWarnings = FALSE, recursive = TRUE)

# Chemin et colonnes de ta table de DEGs RNA-seq.
# RNA_FILE       <- file.path(PROJECT_ROOT, "rnaseq", "DEGs.csv")  # <-- À MODIFIER
RNA_FILE       <- "/home/apuhart/Downloads/Bulk-RNAseq/results_test/fibro/control_vs_FA/fibro_control_vs_FA_DESeq2_results.csv"
RNA_COL_SYMBOL <- "GENENAME"   # nom de la colonne symbole dans ton CSV
RNA_COL_LFC    <- "log2FoldChange"
RNA_COL_PADJ   <- "padj"

# Sens attendu par marque
MARK_TYPE <- c(H3K9me3 = "repressive", H3K27me3 = "repressive",
               H3K27ac = "active",     PolII    = "active")

# Seuils de significativité pour la classification
LFC_TH  <- 1       # |log2FC| minimal (ChIP et RNA)
FDR_TH  <- 0.05    # FDR ChIP
PADJ_TH <- 0.05    # padj RNA

# ---------------------------------------------------------------------------
# LECTURE RNA-seq
# ---------------------------------------------------------------------------
if (!file.exists(RNA_FILE))
  stop(sprintf("Fichier RNA-seq introuvable : %s (modifie RNA_FILE).", RNA_FILE))

rna_raw <- read_csv(RNA_FILE, show_col_types = FALSE)
miss <- setdiff(c(RNA_COL_SYMBOL, RNA_COL_LFC, RNA_COL_PADJ), names(rna_raw))
if (length(miss) > 0)
  stop("Colonnes RNA-seq manquantes : ", paste(miss, collapse = ", "))

rna <- rna_raw %>%
  transmute(GeneSymbol = .data[[RNA_COL_SYMBOL]],
            logFC_RNA  = as.numeric(.data[[RNA_COL_LFC]]),
            padj_RNA   = as.numeric(.data[[RNA_COL_PADJ]])) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "")

# ---------------------------------------------------------------------------
# FUSION POUR CHAQUE TABLE GENE-LEVEL PRODUITE PAR LE SCRIPT 04
# ---------------------------------------------------------------------------
gene_files <- list.files(ANN_DIR, pattern = "_gene_level\\.csv$", full.names = TRUE)
if (length(gene_files) == 0)
  stop("Aucune table *_gene_level.csv dans ", ANN_DIR, " — lance d'abord 04.")

classify <- function(df, mark) {
  type <- MARK_TYPE[[mark]]
  if (is.null(type)) type <- "active"   # défaut prudent
  sig_chip <- !is.na(df$FDR_chip) & df$FDR_chip < FDR_TH & abs(df$logFC_chip) >= LFC_TH
  sig_rna  <- !is.na(df$padj_RNA) & df$padj_RNA < PADJ_TH & abs(df$logFC_RNA) >= LFC_TH

  if (type == "repressive") {
    df$category <- dplyr::case_when(
      sig_chip & sig_rna & df$logFC_chip > 0 & df$logFC_RNA < 0 ~ "gain_mark_down_expr",
      sig_chip & sig_rna & df$logFC_chip < 0 & df$logFC_RNA > 0 ~ "loss_mark_up_expr",
      sig_chip & sig_rna                                        ~ "concordant_unexpected",
      TRUE                                                      ~ "ns_or_other")
  } else {
    df$category <- dplyr::case_when(
      sig_chip & sig_rna & df$logFC_chip > 0 & df$logFC_RNA > 0 ~ "gain_mark_up_expr",
      sig_chip & sig_rna & df$logFC_chip < 0 & df$logFC_RNA < 0 ~ "loss_mark_down_expr",
      sig_chip & sig_rna                                        ~ "discordant_unexpected",
      TRUE                                                      ~ "ns_or_other")
  }
  df
}

for (f in gene_files) {
  base   <- str_remove(basename(f), "_gene_level\\.csv$")  # <MARK>__<CONTRASTE>
  mark   <- str_split(base, "__", simplify = TRUE)[1]
  contrast <- str_split(base, "__", simplify = TRUE)[2]

  chip <- read_csv(f, show_col_types = FALSE) %>%
    transmute(GeneID, GeneSymbol,
              signal_chip = mean_signal,
              logFC_chip  = logFC,
              FDR_chip    = FDR)

  merged <- rna %>%
    inner_join(chip, by = "GeneSymbol") %>%   # gènes communs RNA & ChIP
    classify(mark)

  # Renommage final avec le nom de la marque (table "cible")
  out <- merged %>%
    rename(!!paste0("signal_", mark) := signal_chip,
           !!paste0("logFC_",  mark) := logFC_chip,
           !!paste0("FDR_",    mark) := FDR_chip) %>%
    arrange(category, .data[[paste0("FDR_", mark)]])

  outfile <- file.path(INT_DIR, sprintf("integration_%s__%s.csv", mark, contrast))
  write_csv(out, outfile)
  message(sprintf("  %s : %d gènes communs -> %s",
                  base, nrow(out), basename(outfile)))

  # Sous-table des candidats prioritaires (catégories attendues)
  cand <- out %>% filter(str_detect(category, "^(gain|loss)_"))
  if (nrow(cand) > 0) {
    write_csv(cand, file.path(INT_DIR,
              sprintf("candidates_%s__%s.csv", mark, contrast)))
    message(sprintf("      -> %d candidats anti/corrélés", nrow(cand)))
  }
}

# ---------------------------------------------------------------------------
# (OPTIONNEL) DÉTECTION DE RESTAURATION DANS LES CELLULES ÉDITÉES
# ---------------------------------------------------------------------------
# Principe : un gène est "restauré" si le contraste maladie (ex: FRDA_vs_Control)
# montre un changement et le contraste d'édition (ex: Edited_vs_Unedited) montre
# le changement INVERSE pour la même marque.
# Renseigne ci-dessous les libellés EXACTS de tes contrastes (visibles dans les
# noms de fichiers integration_<MARK>__<CONTRASTE>.csv) puis décommente.
#
# DISEASE_CONTRAST <- "FRDA_vs_Control"
# EDIT_CONTRAST    <- "Edited_vs_Unedited"
#
# for (mark in names(MARK_TYPE)) {
#   fd <- file.path(INT_DIR, sprintf("integration_%s__%s.csv", mark, DISEASE_CONTRAST))
#   fe <- file.path(INT_DIR, sprintf("integration_%s__%s.csv", mark, EDIT_CONTRAST))
#   if (!file.exists(fd) || !file.exists(fe)) next
#   d <- read_csv(fd, show_col_types = FALSE); e <- read_csv(fe, show_col_types = FALSE)
#   lfc <- paste0("logFC_", mark)
#   j <- inner_join(d, e, by = "GeneSymbol", suffix = c("_disease", "_edit"))
#   j$restored <- sign(j[[paste0(lfc, "_disease")]]) != sign(j[[paste0(lfc, "_edit")]]) &
#                 sign(j$logFC_RNA_disease) != sign(j$logFC_RNA_edit)
#   write_csv(filter(j, restored),
#             file.path(INT_DIR, sprintf("restoration_%s.csv", mark)))
# }

message("\n==> [05] Terminé. Tables d'intégration RNA/ChIP écrites dans results/integration/")
