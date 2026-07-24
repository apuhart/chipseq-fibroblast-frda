#!/usr/bin/env Rscript
# =============================================================================
# 13_pathway_enrichment.R — Analyse de voies (GO / KEGG) sur les listes ChIP
# -----------------------------------------------------------------------------
# Pour chaque marque, définit une liste de gènes "candidats" selon le sens
# attendu chez les FRDA, puis teste l'enrichissement fonctionnel (GO Biological
# Process + KEGG). Agrège le signal -> plus robuste que l'analyse gène-par-gène
# quand la puissance individuelle est faible (n=3).
#
#   marque RÉPRESSIVE (H3K9me3, H3K27me3) : candidats = GAIN  (logFC_chip > TH)
#   marque ACTIVE (H3K27ac, PolII)        : candidats = PERTE (logFC_chip < -TH)
#
# Sorties : results/pathways/<MARK>_GO_BP.csv , <MARK>_KEGG.csv , barplots PDF
# EXPLORATOIRE : seuils volontairement souples (TH=0.5), à interpréter avec prudence.
# =============================================================================

## Installe clusterProfiler si absent (peut prendre du temps la 1re fois)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(dplyr); library(readr)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
INT <- file.path(PROJECT_ROOT, "results", "integration")
OUT <- file.path(PROJECT_ROOT, "results", "pathways")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MARK_TYPE <- c(H3K9me3="repressive", H3K27me3="repressive",
               H3K27ac="active", PolII="active")
TH <- 0.5   # seuil souple sur logFC ChIP (exploratoire)

files <- list.files(INT, pattern = "^integration_.*\\.csv$", full.names = TRUE)
files <- files[!grepl("FXN", files)]

top_barplot <- function(df, title, file) {
  if (is.null(df) || nrow(df) == 0) return(invisible())
  df <- head(df[order(df$p.adjust), ], 15)
  pdf(file, width = 9, height = 6); par(mar = c(5, 22, 4, 2))
  barplot(rev(-log10(df$p.adjust)), horiz = TRUE, names.arg = rev(df$Description),
          las = 1, cex.names = 0.7, col = "#4C72B0",
          xlab = "-log10(p.adjust)", main = title)
  dev.off()
}

for (f in files) {
  mark <- sub("^integration_", "", sub("\\.csv$", "", basename(f)))
  d <- read.csv(f, stringsAsFactors = FALSE)
  lfc_col <- paste0("logFC_", mark)
  if (!lfc_col %in% names(d)) next
  d$lfc <- suppressWarnings(as.numeric(d[[lfc_col]]))
  d$ENTREZID <- as.character(d$ENTREZID)

  type <- MARK_TYPE[[mark]]; if (is.null(type)) type <- "active"
  fg <- if (type == "repressive") d$ENTREZID[d$lfc >  TH & !is.na(d$lfc)]
        else                       d$ENTREZID[d$lfc < -TH & !is.na(d$lfc)]
  fg <- unique(na.omit(fg))
  universe <- unique(na.omit(d$ENTREZID))
  cat(sprintf("\n[%s | %s] %d gènes candidats / %d univers\n",
              mark, type, length(fg), length(universe)))
  if (length(fg) < 10) { cat("  trop peu de candidats, marque ignorée.\n"); next }

  go <- tryCatch(enrichGO(fg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                          ont = "BP", universe = universe,
                          pAdjustMethod = "BH", pvalueCutoff = 0.1,
                          qvalueCutoff = 0.2, readable = TRUE),
                 error = function(e) NULL)
  if (!is.null(go) && nrow(as.data.frame(go))) {
    write_csv(as.data.frame(go), file.path(OUT, paste0(mark, "_GO_BP.csv")))
    top_barplot(as.data.frame(go), paste("GO BP —", mark),
                file.path(OUT, paste0(mark, "_GO_BP.pdf")))
    cat("  GO BP : ", nrow(as.data.frame(go)), " termes\n")
  } else cat("  GO BP : aucun terme significatif\n")

  kg <- tryCatch(enrichKEGG(fg, organism = "hsa", universe = universe,
                            pvalueCutoff = 0.1),
                 error = function(e) NULL)
  if (!is.null(kg) && nrow(as.data.frame(kg))) {
    write_csv(as.data.frame(kg), file.path(OUT, paste0(mark, "_KEGG.csv")))
    cat("  KEGG : ", nrow(as.data.frame(kg)), " voies\n")
  } else cat("  KEGG : aucune voie significative\n")
}
cat("\n==> [13] Terminé. Résultats dans results/pathways/\n")
cat("Regarde en priorité : métabolisme mitochondrial, clusters Fe-S, stress oxydatif.\n")
