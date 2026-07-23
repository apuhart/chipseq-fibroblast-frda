#!/usr/bin/env Rscript
# =============================================================================
# 07_fxn_analysis.R — Analyse quantitative ciblée du locus FXN
# -----------------------------------------------------------------------------
# Lit results/fxn/fxn_counts.csv (produit par 06), calcule :
#   - CPM ChIP et input
#   - enrichissement log2((ChIP_CPM+1)/(input_CPM+1))  [signal corrigé du bruit]
#   - moyenne FRDA vs Control par marque et par sous-région
#   - test t (n=3 vs 3 : indicatif, faible puissance)
# Produit :
#   results/fxn/fxn_enrichment.csv          (toutes les valeurs)
#   results/fxn/fxn_summary_CvsF.csv        (résumé C vs F + p-value)
#   results/fxn/fxn_barplots.pdf            (figures)
# Base R uniquement (pas de dépendance ggplot2).
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
FXN_DIR <- file.path(PROJECT_ROOT, "results", "fxn")
META    <- file.path(PROJECT_ROOT, "metadata", "samples_metadata.csv")

counts <- read.csv(file.path(FXN_DIR, "fxn_counts.csv"), stringsAsFactors = FALSE)
meta   <- read.csv(META, stringsAsFactors = FALSE)

# Condition par échantillon
counts$Condition <- meta$Condition[match(counts$sample, meta$sample)]

# CPM et enrichissement
num <- function(x) suppressWarnings(as.numeric(x))
counts$chip_cpm  <- num(counts$chip_count)  / num(counts$chip_libsize)  * 1e6
counts$input_cpm <- num(counts$input_count) / num(counts$input_libsize) * 1e6
counts$enrichment_log2 <- log2((counts$chip_cpm + 1) / (counts$input_cpm + 1))

write.csv(counts, file.path(FXN_DIR, "fxn_enrichment.csv"), row.names = FALSE)

# Résumé FRDA vs Control par marque x région
marks   <- unique(counts$mark)
regions <- unique(counts$region)
res <- data.frame()
for (mk in marks) for (rg in regions) {
  sub <- counts[counts$mark == mk & counts$region == rg, ]
  ctrl <- sub$enrichment_log2[sub$Condition == "Control"]
  frda <- sub$enrichment_log2[sub$Condition == "FRDA"]
  if (length(ctrl) < 2 || length(frda) < 2) next
  p <- tryCatch(t.test(frda, ctrl)$p.value, error = function(e) NA)
  res <- rbind(res, data.frame(
    mark = mk, region = rg,
    mean_Control = round(mean(ctrl), 3),
    mean_FRDA    = round(mean(frda), 3),
    diff_FRDA_minus_Control = round(mean(frda) - mean(ctrl), 3),
    p_value = signif(p, 3)
  ))
}
write.csv(res, file.path(FXN_DIR, "fxn_summary_CvsF.csv"), row.names = FALSE)
cat("\n===== Résumé FXN : enrichissement (log2 ChIP/input), FRDA vs Control =====\n")
print(res, row.names = FALSE)

# --- Figures : un panneau par marque, barres = sous-régions, C vs F ---
pdf(file.path(FXN_DIR, "fxn_barplots.pdf"), width = 9, height = 7)
op <- par(mfrow = c(2, 2), mar = c(7, 4, 3, 1))
cols <- c(Control = "#4C72B0", FRDA = "#C44E52")
for (mk in marks) {
  sub <- counts[counts$mark == mk, ]
  # matrice région x condition (moyenne)
  m <- tapply(sub$enrichment_log2,
              list(sub$region, sub$Condition), mean)
  m <- m[, c("Control", "FRDA"), drop = FALSE]
  bp <- barplot(t(m), beside = TRUE, col = cols, las = 2,
                ylab = "log2(ChIP/input)", main = paste("FXN —", mk),
                legend.text = TRUE,
                args.legend = list(x = "topright", bty = "n"))
  # points individuels superposés (2 colonnes par région : Control, FRDA)
  for (i in seq_len(nrow(m))) {
    rg <- rownames(m)[i]
    for (k in seq_along(c("Control", "FRDA"))) {
      cond <- c("Control", "FRDA")[k]
      vals <- sub$enrichment_log2[sub$region == rg & sub$Condition == cond]
      if (length(vals)) points(rep(bp[k, i], length(vals)), vals, pch = 19, cex = 0.8)
    }
  }
}
par(op); dev.off()
cat("\nFigures -> ", file.path(FXN_DIR, "fxn_barplots.pdf"), "\n")
cat("Interprétation attendue (modèle FRDA) :\n",
    " H3K9me3  intron1/promoteur : FRDA > Control (gain hétérochromatine)\n",
    " H3K27ac / PolII promoteur  : FRDA < Control (perte marque active)\n")
