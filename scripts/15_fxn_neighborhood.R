#!/usr/bin/env Rscript
# =============================================================================
# 15_fxn_neighborhood.R — Propagation de l'hétérochromatine autour de FXN
# -----------------------------------------------------------------------------
# L'expansion (GAA)n peut nucléer du H3K9me3 qui se PROPAGE le long du
# chromosome. Ce script regarde, dans une fenêtre autour de FXN (chr9), si le
# gain de H3K9me3 (et la perte de H3K27ac) chez les FRDA décroît avec la
# distance au locus -> signature de propagation.
#
# Utilise les tables *_regions.csv (qui contiennent les coordonnées).
# Sorties : results/figures/FXN_neighborhood.pdf
#           results/figures/FXN_neighborhood_genes.csv
# =============================================================================
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
DIFF <- file.path(PROJECT_ROOT, "results", "diffbind")
OUT  <- file.path(PROJECT_ROOT, "results", "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FXN_CHR <- "chr9"; FXN_TSS <- 69035563
WINDOW  <- 3e6   # +/- 3 Mb autour de FXN

load_regions <- function(mark) {
  f <- file.path(DIFF, paste0(mark, "__FRDA-vs-Control_regions.csv"))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d <- d[d$seqnames == FXN_CHR, ]
  mid <- (d$start + d$end) / 2
  d$dist_kb <- (mid - FXN_TSS) / 1e3
  d <- d[abs(d$dist_kb) <= WINDOW / 1e3, ]
  d$logFC <- suppressWarnings(as.numeric(d$logFC))
  d[is.finite(d$logFC), ]
}

marks <- c("H3K9me3", "H3K27ac", "H3K27me3", "PolII")
pdf(file.path(OUT, "FXN_neighborhood.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

all_genes <- data.frame()
for (mk in marks) {
  d <- load_regions(mk)
  if (is.null(d) || nrow(d) == 0) {
    plot.new(); title(paste(mk, "- aucune région chr9 dans la fenêtre")); next
  }
  plot(d$dist_kb, d$logFC, pch = 16, cex = 0.6,
       col = ifelse(d$logFC > 0, "#C44E52", "#4C72B0"),
       xlab = "Distance à FXN (kb)", ylab = paste0("logFC ", mk),
       main = paste0(mk, " autour de FXN (+/-3 Mb)"))
  abline(h = 0, lty = 3); abline(v = 0, lty = 2, col = "blue")
  text(0, max(d$logFC, na.rm = TRUE), "FXN", col = "blue", pos = 4, cex = 0.8)
  # tendance lissée
  if (nrow(d) > 5) {
    lo <- tryCatch(loess(logFC ~ dist_kb, data = d, span = 0.6),
                   error = function(e) NULL)
    if (!is.null(lo)) {
      xx <- seq(min(d$dist_kb), max(d$dist_kb), length.out = 100)
      lines(xx, predict(lo, xx), lwd = 2)
    }
  }
  keep <- !is.na(d$SYMBOL) & d$SYMBOL != ""
  if (any(keep)) all_genes <- rbind(all_genes, data.frame(
    mark = mk, gene = d$SYMBOL[keep], dist_kb = round(d$dist_kb[keep]),
    logFC = round(d$logFC[keep], 3), FDR = signif(as.numeric(d$FDR[keep]), 3)))
}
dev.off()

if (nrow(all_genes)) {
  all_genes <- all_genes[order(all_genes$mark, all_genes$dist_kb), ]
  write.csv(all_genes, file.path(OUT, "FXN_neighborhood_genes.csv"), row.names = FALSE)
}
cat("Figure -> ", file.path(OUT, "FXN_neighborhood.pdf"), "\n")
cat("Gènes du voisinage -> ", file.path(OUT, "FXN_neighborhood_genes.csv"), "\n")
cat("\nLecture : si le gain de H3K9me3 (points rouges) est concentré PRÈS de\n",
    "FXN (distance ~0) et décroît avec la distance, c'est compatible avec une\n",
    "propagation locale de l'hétérochromatine.\n")
