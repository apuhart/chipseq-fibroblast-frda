#!/usr/bin/env Rscript
# =============================================================================
# 11_integration_plots.R — Visualisation du croisement RNA-seq <-> ChIP-seq
# -----------------------------------------------------------------------------
# Pour chaque marque, lit results/integration/integration_<MARK>.csv et trace
# un nuage de points : x = logFC_RNA, y = logFC_<MARK>.
#   - Chaque point = un gène.
#   - Les quadrants indiquent la relation chromatine/expression.
#   - FXN est repéré et étiqueté.
#   - Corrélation de Spearman affichée dans le titre.
# Couleurs des quadrants selon le sens attendu :
#   marque RÉPRESSIVE (H3K9me3, H3K27me3) : anticorrélation attendue
#       -> quadrants (logFC_chip>0 & RNA<0) et (chip<0 & RNA>0) = "attendu"
#   marque ACTIVE (H3K27ac, PolII) : corrélation positive attendue
#       -> quadrants (chip>0 & RNA>0) et (chip<0 & RNA<0) = "attendu"
# Sorties : results/figures/integration_scatter.pdf
#           results/figures/integration_quadrant_counts.csv
# Base R uniquement.
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
INT <- file.path(PROJECT_ROOT, "results", "integration")
OUT <- file.path(PROJECT_ROOT, "results", "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MARK_TYPE <- c(H3K9me3="repressive", H3K27me3="repressive",
               H3K27ac="active", PolII="active")

files <- list.files(INT, pattern = "^integration_.*\\.csv$", full.names = TRUE)
files <- files[!grepl("FXN", files)]
if (length(files) == 0) stop("Aucun fichier integration_*.csv. Lance 05b d'abord.")

pdf(file.path(OUT, "integration_scatter.pdf"), width = 9, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

counts_tab <- data.frame()

for (f in files) {
  mark  <- sub("^integration_", "", sub("\\.csv$", "", basename(f)))
  d <- read.csv(f, stringsAsFactors = FALSE)
  lfc_col <- paste0("logFC_", mark)
  if (!all(c("logFC_RNA", lfc_col) %in% names(d))) next
  
  x <- suppressWarnings(as.numeric(d$logFC_RNA))   # RNA
  y <- suppressWarnings(as.numeric(d[[lfc_col]]))   # ChIP
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]; sym <- d$GeneSymbol[keep]
  if (length(x) < 3) next
  
  type <- MARK_TYPE[[mark]]; if (is.null(type)) type <- "active"
  
  # Point "attendu" ou non selon le type de marque
  if (type == "repressive") {
    expected <- (y > 0 & x < 0) | (y < 0 & x > 0)   # anticorrélation
  } else {
    expected <- (y > 0 & x > 0) | (y < 0 & x < 0)   # corrélation positive
  }
  col <- ifelse(expected, rgb(0.77, 0.31, 0.32, 0.5), rgb(0.5, 0.5, 0.5, 0.3))
  
  rho <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
  
  plot(x, y, pch = 16, cex = 0.4, col = col,
       xlab = "logFC RNA (FRDA vs Control)",
       ylab = paste0("logFC ", mark),
       main = sprintf("%s — %s\nSpearman rho = %.3f (n=%d)",
                      mark, type, rho, length(x)))
  abline(h = 0, v = 0, lty = 3, col = "grey40")
  
  # FXN
  i <- which(sym == "FXN")
  if (length(i)) {
    points(x[i], y[i], pch = 1, cex = 2.2, lwd = 2, col = "blue")
    text(x[i], y[i], "FXN", pos = 3, col = "blue", cex = 0.95)
  }
  
  counts_tab <- rbind(counts_tab, data.frame(
    mark = mark, type = type, n_genes = length(x),
    n_expected_quadrant = sum(expected), spearman_rho = round(rho, 3)))
}
dev.off()

write.csv(counts_tab, file.path(OUT, "integration_quadrant_counts.csv"),
          row.names = FALSE)
cat("Scatter -> ", file.path(OUT, "integration_scatter.pdf"), "\n\n")
cat("Résumé par marque :\n"); print(counts_tab, row.names = FALSE)
cat("\nLecture : pour une marque répressive, rho NÉGATIF = anticorrélation\n",
    "globale chromatine/expression (cohérent avec le modèle).\n",
    "Avec n=3 réplicats, attends-toi à des rho faibles : regarde la TENDANCE,\n",
    "pas une corrélation forte.\n")