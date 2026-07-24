#!/usr/bin/env Rscript
# =============================================================================
# 04_diffbind_chipseeker.R
# -----------------------------------------------------------------------------
# Pour chaque marque :
#   1. Construit une feuille d'échantillons DiffBind (manifest + métadonnées).
#   2. Compte les lectures sur le consensus (dba.count) — summits recentrés
#      pour les marques narrow, conservés tels quels pour les broad.
#   3. Analyse différentielle DESeq2 (dba.contrast / dba.analyze) pour chaque
#      paire de conditions ayant >=2 réplicats.
#   4. Annote les régions différentielles avec ChIPseeker (promoteur, exon,
#      intron, intergénique, gène le plus proche).
#   5. Agrège au niveau gène : GeneID, GeneSymbol, signal moyen, logFC, FDR.
#
# Sorties (toutes en CSV) :
#   results/diffbind/<MARK>__<CONTRASTE>_regions.csv     (régions + annotation)
#   results/annotation/<MARK>__<CONTRASTE>_gene_level.csv
#   results/annotation/<MARK>_consensus_annotation.csv   (descriptif)
#
# Lancer depuis la racine du projet :  Rscript scripts/04_diffbind_chipseeker.R
# =============================================================================

suppressPackageStartupMessages({
  library(DiffBind)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)   # <-- adapter si autre build
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(dplyr)
  library(readr)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
META_DIR <- file.path(PROJECT_ROOT, "metadata")
RES      <- file.path(PROJECT_ROOT, "results")
DIFF_DIR <- file.path(RES, "diffbind")
ANN_DIR  <- file.path(RES, "annotation")
CONS_DIR <- file.path(RES, "consensus")
dir.create(DIFF_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(ANN_DIR,  showWarnings = FALSE, recursive = TRUE)

CONTRAST_FACTOR <- "Condition"      # colonne des métadonnées définissant les groupes
TSS_REGION      <- c(-3000, 3000)   # définition du promoteur pour ChIPseeker
TXDB            <- TxDb.Hsapiens.UCSC.hg38.knownGene
MIN_REPLICATES  <- 2                # minimum requis par groupe pour le test

# Choix de l'agrégation pic -> gène (documenté) :
#   "best_fdr" : on garde, par gène, la région la plus significative (FDR min).
#   "max_abs"  : on garde la région au |logFC| le plus fort.
AGGREGATION <- "best_fdr"

# ---------------------------------------------------------------------------
# LECTURE DES ENTREES
# ---------------------------------------------------------------------------
manifest <- read_csv(file.path(META_DIR, "chip_manifest_with_peaks.csv"),
                     show_col_types = FALSE)
meta     <- read_csv(file.path(META_DIR, "samples_metadata.csv"),
                     show_col_types = FALSE)

if (!CONTRAST_FACTOR %in% names(meta))
  stop(sprintf("Colonne '%s' absente de samples_metadata.csv", CONTRAST_FACTOR))

ss <- manifest %>% left_join(meta, by = "sample")

# Garde-fou : conditions renseignées ?
if (all(is.na(ss[[CONTRAST_FACTOR]]) | ss[[CONTRAST_FACTOR]] == "")) {
  stop("La colonne Condition est vide dans samples_metadata.csv. ",
       "Complète-la (ex: FRDA / Control / Edited) avant de relancer.")
}

# ---------------------------------------------------------------------------
# FONCTIONS UTILITAIRES
# ---------------------------------------------------------------------------

# Construit la feuille d'échantillons DiffBind pour un sous-ensemble (1 marque).
make_sheet <- function(df) {
  data.frame(
    SampleID   = paste(df$sample, df$mark, df$Replicate, sep = "_"),
    Tissue     = "Fibroblast",
    Factor     = df$mark,
    Condition  = df[[CONTRAST_FACTOR]],
    Treatment  = ifelse(is.na(df$Treatment), "NA", df$Treatment),
    Replicate  = df$Replicate,
    bamReads   = df$chip_bam,
    ControlID  = paste0(df$sample, "_input"),
    bamControl = df$input_bam,
    Peaks      = df$peaks,
    # narrowPeak -> "narrow" ; broadPeak (BED6+3) lu comme "bed"
    PeakCaller = ifelse(df$peaktype == "narrow", "narrow", "bed"),
    stringsAsFactors = FALSE
  )
}

# Simplifie la catégorie d'annotation ChIPseeker.
simplify_anno <- function(x) {
  dplyr::case_when(
    grepl("^Promoter",      x) ~ "Promoter",
    grepl("^Exon",          x) ~ "Exon",
    grepl("^Intron",        x) ~ "Intron",
    grepl("Intergenic",     x) ~ "Intergenic",
    grepl("UTR",            x) ~ "UTR",
    grepl("^Downstream",    x) ~ "Downstream",
    TRUE                       ~ "Other"
  )
}

# Annote un GRanges avec ChIPseeker et renvoie un data.frame propre.
annotate_gr <- function(gr) {
  if (length(gr) == 0) return(data.frame())
  # Harmonise le style de chromosomes vers UCSC (chr1, chr2, ...) pour matcher le TxDb.
  tryCatch(seqlevelsStyle(gr) <- "UCSC", error = function(e) {})
  ap <- annotatePeak(gr, TxDb = TXDB, tssRegion = TSS_REGION,
                     annoDb = "org.Hs.eg.db", verbose = FALSE)
  df <- as.data.frame(ap)
  df$annotation_simple <- simplify_anno(df$annotation)
  df
}

# Agrège les régions annotées au niveau gène.
aggregate_genes <- function(df) {
  needed <- c("geneId", "SYMBOL", "Conc", "logFC", "FDR")
  if (!all(needed %in% names(df))) return(data.frame())
  df %>%
    filter(!is.na(SYMBOL), SYMBOL != "") %>%
    group_by(GeneID = geneId, GeneSymbol = SYMBOL) %>%
    summarise(
      n_regions   = dplyr::n(),
      mean_signal = mean(Conc, na.rm = TRUE),
      logFC = if (AGGREGATION == "max_abs")
                logFC[which.max(abs(logFC))] else logFC[which.min(FDR)],
      FDR   = if (AGGREGATION == "max_abs")
                FDR[which.max(abs(logFC))]   else min(FDR, na.rm = TRUE),
      .groups = "drop"
    )
}

# ---------------------------------------------------------------------------
# BOUCLE PRINCIPALE PAR MARQUE
# ---------------------------------------------------------------------------
marks <- unique(ss$mark)
message(sprintf("Marques détectées : %s", paste(marks, collapse = ", ")))

for (mk in marks) {
  message(sprintf("\n================ %s ================", mk))
  sub   <- ss %>% filter(mark == mk)
  sheet <- make_sheet(sub)

  # ---- Annotation descriptive du consensus (toujours produite) ----
  cons_bed <- file.path(CONS_DIR, paste0(mk, "_consensus_union.bed"))
  if (file.exists(cons_bed) && file.info(cons_bed)$size > 0) {
    cons <- read_tsv(cons_bed, col_names = c("chr", "start", "end"),
                     show_col_types = FALSE)
    cons_gr <- GRanges(cons$chr, IRanges(cons$start + 1, cons$end))  # BED -> 1-based
    cdf <- annotate_gr(cons_gr)
    if (nrow(cdf) > 0)
      write_csv(cdf, file.path(ANN_DIR, paste0(mk, "_consensus_annotation.csv")))
  }

  # ---- Vérification des conditions / réplicats ----
  ctab <- table(sheet$Condition)
  if (length(ctab) < 2 || any(ctab < MIN_REPLICATES)) {
    message(sprintf("  [SKIP différentiel] %s : besoin de >=%d réplicats dans >=2 conditions. Répartition : %s",
                    mk, MIN_REPLICATES,
                    paste(names(ctab), ctab, sep = "=", collapse = ", ")))
    next
  }

  # ---- DiffBind ----
  is_broad <- unique(sub$peaktype) == "broad"
  summ <- if (isTRUE(is_broad)) FALSE else 200   # ne pas recentrer les broad

  dbo <- dba(sampleSheet = sheet)
  # bUseSummarizeOverlaps=TRUE requiert des BAM indexés (.bai).
  dbo <- dba.count(dbo, bUseSummarizeOverlaps = TRUE, summits = summ)
  dbo <- dba.normalize(dbo)                       # normalisation par défaut (lib size)
  dbo <- dba.contrast(dbo, categories = DBA_CONDITION, minMembers = MIN_REPLICATES)
  dbo <- dba.analyze(dbo, method = DBA_DESEQ2)

  contrasts_info <- dba.show(dbo, bContrasts = TRUE)
  if (is.null(contrasts_info) || nrow(contrasts_info) == 0) {
    message(sprintf("  [%s] Aucun contraste exploitable.", mk)); next
  }

  # ---- Boucle sur chaque contraste ----
  for (ci in seq_len(nrow(contrasts_info))) {
    g1 <- as.character(contrasts_info$Group[ci])
    g2 <- as.character(contrasts_info$Group2[ci])
    clabel <- gsub("[^A-Za-z0-9]+", "-", paste0(g1, "_vs_", g2))

    # th=1 -> récupère TOUTES les régions (l'utilisateur filtre ensuite).
    rep_gr <- tryCatch(dba.report(dbo, contrast = ci, method = DBA_DESEQ2, th = 1),
                       error = function(e) { message("    report KO: ", e$message); NULL })
    if (is.null(rep_gr) || length(rep_gr) == 0) {
      message(sprintf("    [%s | %s] rapport vide.", mk, clabel)); next
    }

    # Uniformise les noms de colonnes (DiffBind nomme 'Fold' et 'p-value').
    mc <- mcols(rep_gr)
    if ("Fold"    %in% names(mc)) names(mc)[names(mc) == "Fold"]    <- "logFC"
    if ("p-value" %in% names(mc)) names(mc)[names(mc) == "p-value"] <- "pvalue"
    mcols(rep_gr) <- mc

    # Annotation ChIPseeker des régions différentielles
    adf <- annotate_gr(rep_gr)
    adf$mark <- mk; adf$contrast <- clabel

    out_regions <- file.path(DIFF_DIR, sprintf("%s__%s_regions.csv", mk, clabel))
    write_csv(adf, out_regions)
    message(sprintf("    [%s | %s] %d régions -> %s",
                    mk, clabel, nrow(adf), basename(out_regions)))

    # Agrégation au niveau gène
    gdf <- aggregate_genes(adf)
    if (nrow(gdf) > 0) {
      gdf$mark <- mk; gdf$contrast <- clabel
      out_genes <- file.path(ANN_DIR, sprintf("%s__%s_gene_level.csv", mk, clabel))
      write_csv(gdf, out_genes)
      message(sprintf("    [%s | %s] %d gènes -> %s",
                      mk, clabel, nrow(gdf), basename(out_genes)))
    }
  }
}

message("\n==> [04] Terminé. Tables différentielles + niveau gène écrites.")
