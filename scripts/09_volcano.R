#!/usr/bin/env Rscript
# =============================================================================
# 09_volcano.R — Volcano plots (logFC vs significativité) par marque
# -----------------------------------------------------------------------------
# Visualise, pour chaque marque, l'ensemble des régions testées :
#   x = logFC (FRDA vs Control), y = -log10(p-value)
# Met en évidence le peu (ou beaucoup) de régions significatives, et où
# tombe FXN. Base R, sortie PDF unique.
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
DIFF <- file.path(PROJECT_ROOT, "results", "diffbind")
OUT  <- file.path(PROJECT_ROOT, "results", "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

files <- list.files(DIFF, pattern = "_regions\\.csv$", full.names = TRUE)
FDR_TH <- 0.05

pdf(file.path(OUT, "volcano_all_marks.pdf"), width = 9, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (f in files) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  mark <- sub("__.*", "", basename(f))
  d$logFC  <- suppressWarnings(as.numeric(d$logFC))
  d$pvalue <- suppressWarnings(as.numeric(d$pvalue))
  d$FDR    <- suppressWarnings(as.numeric(d$FDR))
  d <- d[is.finite(d$logFC) & is.finite(d$pvalue), ]

  sig <- !is.na(d$FDR) & d$FDR < FDR_TH
  col <- ifelse(sig, "#C44E52", "grey70")

  plot(d$logFC, -log10(d$pvalue), pch = 16, cex = 0.4, col = col,
       xlab = "logFC (FRDA vs Control)", ylab = "-log10(p-value)",
       main = sprintf("%s  (%d régions, %d signif. FDR<%.2f)",
                      mark, nrow(d), sum(sig, na.rm = TRUE), FDR_TH))
  abline(v = 0, lty = 3, col = "grey40")

  # Position de FXN sur ce volcano
  fx <- d[!is.na(d$SYMBOL) & d$SYMBOL == "FXN", ]
  if (nrow(fx)) {
    points(fx$logFC, -log10(fx$pvalue), pch = 1, cex = 2, lwd = 2, col = "blue")
    text(fx$logFC, -log10(fx$pvalue), labels = "FXN",
         pos = 3, col = "blue", cex = 0.9)
  }
}
dev.off()
cat("Volcano -> ", file.path(OUT, "volcano_all_marks.pdf"), "\n")
