#!/usr/bin/env Rscript
# =============================================================================
# 04b_recompute_gene_level.R — CORRECTIF du logFC ChIP + tables gene-level
# -----------------------------------------------------------------------------
# PROBLÈME CORRIGÉ ICI
#   DiffBind 3.x renvoie par défaut un `Fold` SHRUNKÉ (apeglm). Avec une
#   dispersion élevée (n=3/groupe), le shrinkage écrase le fold à ~0 :
#       Conc_FRDA = 5.18 , Conc_Control = 3.58  -> vrai fold ≈ 1.60
#       Fold (shrunké)                          -> 0.0078
#   Conséquence : max|logFC| ≈ 0.01 génome-entier, aucun filtre sur le fold
#   ne peut passer, et TOUTE analyse basée sur l'amplitude est fausse.
#
#   Correctif : logFC <- Conc_<condB> - Conc_<condA>  (MLE non shrunké).
#   Les colonnes Conc_* sont déjà dans les _regions.csv -> AUCUN besoin de
#   relancer DiffBind.
#
# CAVEAT : le fold non shrunké est bruité pour les régions à faible couverture
#   (c'est précisément ce que le shrinkage voulait corriger). Il est donc juste
#   pour : les analyses de SIGNE, les tendances moyennées, les stratifications.
#   Il ne l'est PAS pour : classer les gènes par |fold| sans tenir compte du FDR.
#
# EN PLUS
#   Classe chaque région (Promoter / Genic / Distal) depuis l'annotation
#   ChIPseeker, et produit des tables gene-level PAR CLASSE + une table "All".
#   Cela permet la stratification promoteur vs distal en aval.
#
# Entrée  : results/diffbind/<MARK>__<CONTRASTE>_regions.csv
# Sortie  : results/annotation/<MARK>__<CONTRASTE>_gene_level_fixed.csv
# Usage   : Rscript 04b_recompute_gene_level.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
DB_DIR  <- file.path(PROJECT_ROOT, "results", "diffbind")
ANN_DIR <- file.path(PROJECT_ROOT, "results", "annotation")
dir.create(ANN_DIR, showWarnings = FALSE, recursive = TRUE)

# Reprend le réglage de 04 : "max_abs" = on garde la région au |logFC| le plus
# fort ; sinon on garde la région au FDR le plus faible.
AGGREGATION <- "max_abs"

# Seuil (bp) définissant le promoteur pour la classification de région
PROMOTER_BP <- 2000

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Classe une annotation ChIPseeker en 3 catégories mécanistiquement distinctes.
# On ne mélange PAS promoteur et enhancer distal : ce sont deux mécanismes.
classify_region <- function(annotation, distanceToTSS = NULL) {
  ann <- as.character(annotation)
  out <- dplyr::case_when(
    str_detect(ann, regex("^Promoter", ignore_case = TRUE))          ~ "Promoter",
    str_detect(ann, regex("5' UTR|3' UTR|Exon|Intron", ignore_case = TRUE)) ~ "Genic",
    str_detect(ann, regex("Distal Intergenic|Downstream", ignore_case = TRUE)) ~ "Distal",
    TRUE ~ "Other")
  # Si distanceToTSS est disponible, il tranche pour le promoteur (plus fiable
  # que la chaîne de caractères, dont le format varie selon les versions).
  if (!is.null(distanceToTSS)) {
    d <- suppressWarnings(as.numeric(distanceToTSS))
    out[!is.na(d) & abs(d) <= PROMOTER_BP] <- "Promoter"
  }
  out
}

# Agrège les régions d'un gène en une ligne.
#
# ATTENTION — PIÈGE dplyr (bug corrigé ici) :
#   summarise() évalue ses expressions SÉQUENTIELLEMENT, et chaque colonne
#   créée MASQUE la variable d'origine pour les expressions suivantes. Écrire
#       logFC      = logFC[which.max(abs(logFC))],   # logFC devient un SCALAIRE
#       FDR        = FDR[which.max(abs(logFC))],     # <- abs(SCALAIRE) !
#       mean_logFC = mean(logFC)                     # <- mean(SCALAIRE) = SCALAIRE
#   produit : un FDR pris sur la 1re région (which.max d'un scalaire = 1),
#   et un mean_logFC identique à logFC. Le test de robustesse renvoyait
#   alors delta = 0 partout, ce qui ressemblait à un succès.
#
#   PARADE : on ne crée JAMAIS de colonne portant le nom d'une variable qu'on
#   doit encore lire. Toutes les statistiques sont calculées sous des noms
#   distincts, et logFC/FDR ne sont assignés qu'à la toute fin, par mutate().
aggregate_gene <- function(d) {
  d %>%
    summarise(
      n_regions     = dplyr::n(),
      mean_signal   = mean(Conc, na.rm = TRUE),
      # --- statistiques sur le VECTEUR de régions (noms distincts) ---
      mean_logFC    = mean(logFC, na.rm = TRUE),
      med_logFC     = median(logFC, na.rm = TRUE),
      min_FDR       = min(FDR, na.rm = TRUE),
      # --- région au |logFC| le plus fort ---
      maxabs_logFC  = logFC[which.max(abs(logFC))],
      maxabs_FDR    = FDR[which.max(abs(logFC))],
      # --- région au meilleur FDR ---
      bestfdr_logFC = logFC[which.min(FDR)],
      bestfdr_FDR   = FDR[which.min(FDR)],
      .groups = "drop") %>%
    # Colonnes canoniques, choisies APRÈS coup selon AGGREGATION.
    mutate(
      logFC = if (AGGREGATION == "max_abs") maxabs_logFC else bestfdr_logFC,
      FDR   = if (AGGREGATION == "max_abs") maxabs_FDR   else bestfdr_FDR)
}

# ---------------------------------------------------------------------------
# Boucle sur les fichiers de régions
# ---------------------------------------------------------------------------
files <- list.files(DB_DIR, pattern = "_regions\\.csv$", full.names = TRUE)
if (!length(files)) stop("Aucun *_regions.csv dans ", DB_DIR)

cat("\n=================== 04b : CORRECTIF logFC ChIP ===================\n\n")

for (f in files) {
  
  tag  <- sub("_regions$", "", tools::file_path_sans_ext(basename(f)))
  mark <- sub("__.*$", "", tag)
  
  d <- read_csv(f, show_col_types = FALSE)
  
  # ---- Repère les colonnes de concentration par groupe ---------------------
  conc_cols <- grep("^Conc_", names(d), value = TRUE)
  if (length(conc_cols) != 2) {
    warning(sprintf("[%s] %d colonnes Conc_* trouvées (2 attendues) : %s — ignoré.",
                    tag, length(conc_cols), paste(conc_cols, collapse = ", ")))
    next
  }
  
  # Le contraste est encodé dans le nom de fichier : <MARK>__<condB>-vs-<condA>
  # DiffBind nomme Conc_<group1> puis Conc_<group2>, et Fold = group1 - group2.
  contrast_part <- sub("^.*__", "", tag)              # ex : "FRDA-vs-Control"
  condB <- sub("-vs-.*$", "", contrast_part)          # ex : "FRDA"
  condA <- sub("^.*-vs-", "", contrast_part)          # ex : "Control"
  
  colB <- paste0("Conc_", condB)
  colA <- paste0("Conc_", condA)
  
  if (!all(c(colB, colA) %in% names(d))) {
    # Repli : on suit l'ordre des colonnes tel que DiffBind les écrit.
    warning(sprintf("[%s] %s / %s introuvables, repli sur l'ordre des colonnes (%s - %s).",
                    tag, colB, colA, conc_cols[1], conc_cols[2]))
    colB <- conc_cols[1]; colA <- conc_cols[2]
  }
  
  # ---- LE CORRECTIF -------------------------------------------------------
  logFC_old <- d$logFC
  d$logFC   <- as.numeric(d[[colB]]) - as.numeric(d[[colA]])
  
  cat(sprintf("[%s]  %s - %s\n", tag, colB, colA))
  cat(sprintf("       max|logFC| :  %.4f (avant, shrunké)  ->  %.3f (après)\n",
              max(abs(logFC_old), na.rm = TRUE),
              max(abs(d$logFC),   na.rm = TRUE)))
  
  # Contrôle de cohérence : le signe doit être préservé par le shrinkage.
  ok <- !is.na(logFC_old) & !is.na(d$logFC) & logFC_old != 0 & d$logFC != 0
  agree <- mean(sign(logFC_old[ok]) == sign(d$logFC[ok]))
  cat(sprintf("       signes conservés vs ancien fold : %.1f%%", 100 * agree))
  cat(if (agree > 0.95) "  (OK : le shrinkage n'inverse pas les signes)\n" else
    "  (!! ATTENTION : incohérence, vérifie l'ordre des groupes)\n")
  
  # ---- Classification des régions -----------------------------------------
  if (!"annotation" %in% names(d)) {
    warning(sprintf("[%s] pas de colonne 'annotation' : classification impossible.", tag))
    d$region_class <- "Other"
  } else {
    d$region_class <- classify_region(
      d$annotation,
      if ("distanceToTSS" %in% names(d)) d$distanceToTSS else NULL)
  }
  
  # ---- Réécrit les régions corrigées --------------------------------------
  d$logFC_shrunk_orig <- logFC_old   # trace de l'ancienne valeur
  write_csv(d, file.path(DB_DIR, paste0(tag, "_regions_fixed.csv")))
  
  # ---- Agrégation gene-level ----------------------------------------------
  needed <- c("geneId", "SYMBOL", "Conc", "logFC", "FDR")
  miss <- setdiff(needed, names(d))
  if (length(miss)) {
    warning(sprintf("[%s] colonnes manquantes : %s — gene-level ignoré.",
                    tag, paste(miss, collapse = ", ")))
    next
  }
  
  base <- d %>% filter(!is.na(geneId), !is.na(SYMBOL))
  
  # Une ligne par gène x classe de région
  by_class <- base %>%
    group_by(GeneID = geneId, GeneSymbol = SYMBOL, region_class) %>%
    aggregate_gene()
  
  # Plus une ligne "All" (toutes régions confondues) = rétro-compatible avec 05b
  all_reg <- base %>%
    group_by(GeneID = geneId, GeneSymbol = SYMBOL) %>%
    aggregate_gene() %>%
    mutate(region_class = "All")
  
  gdf <- bind_rows(all_reg, by_class) %>%
    relocate(GeneID, GeneSymbol, region_class) %>%
    arrange(GeneSymbol, region_class)
  
  out <- file.path(ANN_DIR, sprintf("%s_gene_level_fixed.csv", tag))
  write_csv(gdf, out)
  
  cat(sprintf("       %d gènes | classes : %s\n",
              n_distinct(gdf$GeneID),
              paste(sprintf("%s=%d", names(table(by_class$region_class)),
                            table(by_class$region_class)), collapse = " ")))
  cat(sprintf("       -> %s\n\n", basename(out)))
}

cat("==> [04b] Terminé.\n")
cat("    Prochaine étape : relance 05b en pointant sur *_gene_level_fixed.csv,\n")
cat("    puis 16_integration_analysis.R\n\n")