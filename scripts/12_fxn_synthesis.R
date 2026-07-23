#!/usr/bin/env Rscript
# =============================================================================
# 12_fxn_synthesis.R — Figure de synthèse FXN (ChIP ciblé + RNA-seq)
# -----------------------------------------------------------------------------
# Combine en UNE figure l'histoire mécaniste de FXN dans les FRDA :
#   - enrichissement ChIP ciblé (module 06/07) : différence FRDA - Control
#     par marque et sous-région (promoteur / intron1-GAA / corps du gène) ;
#   - expression RNA-seq (logFC FRDA vs Control) en référence.
# Lecture attendue (modèle Friedreich) :
#   marques actives (H3K27ac, PolII) : barres NÉGATIVES (perte) au promoteur
#   marques répressives (H3K9me3, H3K27me3) : barres POSITIVES (gain) sur GAA
#   expression : fortement négative
# Sortie : results/figures/FXN_synthesis.pdf
# Base R uniquement.
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
FXN_DIR <- file.path(PROJECT_ROOT, "results", "fxn")
OUT     <- file.path(PROJECT_ROOT, "results", "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# --- ChIP ciblé : résumé FRDA vs Control ---
sumf <- file.path(FXN_DIR, "fxn_summary_CvsF.csv")
if (!file.exists(sumf)) stop("Lance d'abord 06 + 07 (fxn_summary_CvsF.csv manquant).")
s <- read.csv(sumf, stringsAsFactors = FALSE)
# colonne d'intérêt : diff_FRDA_minus_Control
s$diff <- as.numeric(s$diff_FRDA_minus_Control)

# --- Expression RNA-seq de FXN (lue depuis l'intégration si dispo) ---
logFC_RNA <- NA
fxnint <- list.files(file.path(PROJECT_ROOT, "results", "integration"),
                     pattern = "integration_.*\\.csv$", full.names = TRUE)
for (f in fxnint) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  r <- d[!is.na(d$GeneSymbol) & d$GeneSymbol == "FXN", ]
  if (nrow(r)) { logFC_RNA <- as.numeric(r$logFC_RNA[1]); break }
}

# Ordonne marques et régions de façon lisible
mark_order   <- c("H3K9me3", "H3K27me3", "H3K27ac", "PolII")
region_order <- c("FXN_promoter", "FXN_intron1_GAA", "FXN_genebody")
region_lab   <- c(FXN_promoter="Promoteur", FXN_intron1_GAA="Intron1/GAA",
                  FXN_genebody="Corps du gène")

# Matrice région x marque (diff FRDA-Control)
m <- matrix(NA, nrow = length(region_order), ncol = length(mark_order),
            dimnames = list(region_lab[region_order], mark_order))
for (i in seq_len(nrow(s))) {
  rg <- s$region[i]; mk <- s$mark[i]
  if (rg %in% region_order && mk %in% mark_order)
    m[region_lab[rg], mk] <- s$diff[i]
}

pdf(file.path(OUT, "FXN_synthesis.pdf"), width = 10, height = 6)
layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))

# --- Panneau 1 : ChIP ciblé ---
par(mar = c(5, 4, 4, 1))
cols_region <- c("#7E57C2", "#C44E52", "#4C72B0")  # promoteur, GAA, corps
bp <- barplot(m, beside = TRUE, col = cols_region, las = 1,
              ylab = "Δ enrichissement log2 (FRDA − Control)",
              main = "FXN — ChIP ciblé (signal corrigé de l'input)",
              ylim = range(c(m, 0), na.rm = TRUE) * 1.3,
              legend.text = rownames(m),
              args.legend = list(x = "topright", bty = "n", cex = 0.9))
abline(h = 0, lwd = 1)
mtext("Répressif (gain attendu chez FRDA)", side = 1, line = 3, at = mean(bp[, 1:2]),
      cex = 0.8, col = "grey30")
mtext("Actif (perte attendue chez FRDA)", side = 1, line = 3, at = mean(bp[, 3:4]),
      cex = 0.8, col = "grey30")

# --- Panneau 2 : expression RNA ---
par(mar = c(5, 4, 4, 2))
barplot(logFC_RNA, col = "#55A868", ylab = "logFC RNA (FRDA vs Control)",
        main = "FXN — Expression", names.arg = "FXN",
        ylim = range(c(logFC_RNA, 0), na.rm = TRUE) * 1.3)
abline(h = 0, lwd = 1)
text(0.7, logFC_RNA, sprintf("%.2f", logFC_RNA),
     pos = if (logFC_RNA < 0) 1 else 3, font = 2)

dev.off()
cat("Figure de synthèse -> ", file.path(OUT, "FXN_synthesis.pdf"), "\n")
cat(sprintf("\nFXN — récapitulatif :\n  Expression RNA logFC = %.2f\n", logFC_RNA))
cat("  ChIP ciblé (Δ FRDA-Control) :\n")
print(round(m, 3))