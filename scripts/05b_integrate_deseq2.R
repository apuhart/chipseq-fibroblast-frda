#!/usr/bin/env Rscript
# =============================================================================
# 05b_integrate_deseq2.R — Intégration RNA-seq (DESeq2) <-> ChIP-seq
# -----------------------------------------------------------------------------
# ADAPTÉ pour un CSV DESeq2 dont la 1re colonne ('gene') contient des SYMBOLES
# (avec quelques identifiants Ensembl en repli quand le symbole est absent).
#
# Étapes :
#   1. Lit le CSV DESeq2 (colonnes log2FoldChange / pvalue / padj).
#   2. Mappe la colonne 'gene' -> Entrez, via SYMBOL *et* ENSEMBL (mixte).
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

# ---- Fichier RNA-seq ----
RNA_FILE <- "/home/apuhart/Downloads/Bulk-RNAseq/results_test/fibro/control_vs_FA/fibro_control_vs_FA_DESeq2_results.csv"

# Sens attendu par marque (répressif : anticorrélation ; actif : corrélation)
MARK_TYPE <- c(H3K9me3="repressive", H3K27me3="repressive",
               H3K27ac="active", PolII="active")

# Seuils pour la colonne 'category'.
# NB : avec n=3/groupe, exiger FDR<0.05 des DEUX côtés est très restrictif.
# Mets SIG_TH plus haut (0.1-0.25) pour un screening ; garde 0.05 pour publier.
LFC_TH <- 1; SIG_TH <- 0.05

# ---------------------------------------------------------------------------
# 1. Lecture RNA-seq (DESeq2)
# ---------------------------------------------------------------------------
rna_raw <- read_csv(RNA_FILE, show_col_types = FALSE)

# La 1re colonne peut s'appeler 'gene', '...1', 'GENEID'… : on la prend par position.
names(rna_raw)[1] <- "GENE_ID"

# Vérifie que les colonnes DESeq2 attendues sont là.
need <- c("log2FoldChange", "pvalue", "padj")
miss <- setdiff(need, names(rna_raw))
if (length(miss))
  stop("Colonnes absentes du CSV RNA : ", paste(miss, collapse = ", "),
       "\nColonnes trouvées : ", paste(names(rna_raw), collapse = ", "))

rna <- rna_raw %>% transmute(
  GENE_ID    = sub("\\..*$", "", as.character(GENE_ID)),  # retire suffixe .version
  logFC_RNA  = as.numeric(log2FoldChange),
  FDR_RNA    = as.numeric(padj),
  PValue_RNA = as.numeric(pvalue)) %>%
  filter(!is.na(GENE_ID), GENE_ID != "")

# ---------------------------------------------------------------------------
# 2. Mapping mixte : SYMBOL et/ou ENSEMBL -> ENTREZID
# ---------------------------------------------------------------------------
# La colonne 'gene' est majoritairement composée de SYMBOLES (APOE, STC2…),
# avec des ENSG en repli quand aucun symbole n'est disponible.
# On mappe donc les deux familles séparément, puis on les recolle.
is_ens <- grepl("^ENSG[0-9]+$", rna$GENE_ID)
cat(sprintf("RNA : %d lignes — %d symboles, %d Ensembl.\n",
            nrow(rna), sum(!is_ens), sum(is_ens)))

safe_select <- function(keys_in, keytype, columns) {
  keys_ok <- intersect(unique(keys_in), AnnotationDbi::keys(org.Hs.eg.db, keytype))
  if (!length(keys_ok)) return(NULL)
  suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = keys_ok,
                          keytype = keytype, columns = columns))
}

# --- famille symboles ---
map_sym <- safe_select(rna$GENE_ID[!is_ens], "SYMBOL", c("ENTREZID", "ENSEMBL"))
if (!is.null(map_sym)) {
  map_sym <- map_sym %>%
    distinct(SYMBOL, .keep_all = TRUE) %>%
    transmute(GENE_ID = SYMBOL, ENTREZID, ENSEMBL, SYMBOL)
} else map_sym <- tibble()

# --- famille Ensembl ---
map_ens <- safe_select(rna$GENE_ID[is_ens], "ENSEMBL", c("ENTREZID", "SYMBOL"))
if (!is.null(map_ens)) {
  map_ens <- map_ens %>%
    distinct(ENSEMBL, .keep_all = TRUE) %>%
    transmute(GENE_ID = ENSEMBL, ENTREZID, ENSEMBL, SYMBOL)
} else map_ens <- tibble()

map_all <- bind_rows(map_sym, map_ens) %>% distinct(GENE_ID, .keep_all = TRUE)
rna <- rna %>% left_join(map_all, by = "GENE_ID")

cat(sprintf("      %d/%d mappés vers Entrez (%.1f%%).\n",
            sum(!is.na(rna$ENTREZID)), nrow(rna),
            100 * mean(!is.na(rna$ENTREZID))))
if (mean(!is.na(rna$ENTREZID)) < 0.5)
  warning("Moins de 50% des gènes mappés — vérifie le contenu de la 1re colonne.")

# ---------------------------------------------------------------------------
# 3-4. Fusion avec chaque table ChIP gene_level
# ---------------------------------------------------------------------------
classify <- function(d, type) {
  sc <- !is.na(d$FDR_chip) & d$FDR_chip < SIG_TH & abs(d$logFC_chip) >= LFC_TH
  sr <- !is.na(d$FDR_RNA)  & d$FDR_RNA  < SIG_TH & abs(d$logFC_RNA)  >= LFC_TH
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
if (!length(files)) stop("Aucune table *_gene_level.csv dans ", ANN_DIR)

for (f in files) {
  mark <- sub("__.*", "", basename(f))
  chip <- read_csv(f, show_col_types = FALSE) %>%
    transmute(ENTREZID = as.character(GeneID), GeneSymbol,
              signal_chip = mean_signal, logFC_chip = logFC, FDR_chip = FDR)
  
  # ---- Garde-fou : logFC ChIP dégénéré ? -----------------------------------
  # DiffBind 3.x renvoie par défaut un Fold *shrunké* qui peut être écrasé
  # à ~0 quand la dispersion est élevée. Dans ce cas aucun gène ne peut
  # passer abs(logFC_chip) >= LFC_TH, et 'category' reste vide.
  max_lfc <- max(abs(chip$logFC_chip), na.rm = TRUE)
  if (is.finite(max_lfc) && max_lfc < 0.1) {
    warning(sprintf(
      "[%s] max|logFC_chip| = %.4g : le fold ChIP est écrasé (shrinkage DiffBind).\n",
      mark, max_lfc),
      "  -> aucun gène ne pourra passer LFC_TH. Recalcule logFC dans 04 via\n",
      "     Conc_<grp1> - Conc_<grp2>, ou utilise dba.report(..., precision=0).",
      call. = FALSE)
  }
  
  merged <- rna %>%
    mutate(ENTREZID = as.character(ENTREZID)) %>%
    filter(!is.na(ENTREZID)) %>%
    inner_join(chip, by = "ENTREZID")
  
  if (nrow(merged) == 0) { cat(sprintf("[%s] aucun gène commun.\n", mark)); next }
  merged$category <- classify(merged, MARK_TYPE[[mark]])
  
  out <- merged %>%
    transmute(GeneSymbol = coalesce(GeneSymbol, SYMBOL, GENE_ID),
              GENE_ID, ENSEMBL, ENTREZID,
              logFC_RNA, PValue_RNA, FDR_RNA,
              !!paste0("signal_", mark) := signal_chip,
              !!paste0("logFC_",  mark) := logFC_chip,
              !!paste0("FDR_",    mark) := FDR_chip,
              category) %>%
    arrange(FDR_RNA)
  
  write_csv(out, file.path(INT_DIR, sprintf("integration_%s.csv", mark)))
  
  n_cat <- sum(out$category != "ns_or_other")
  cat(sprintf("[%s] %d gènes communs, %d catégorisés -> integration_%s.csv\n",
              mark, nrow(out), n_cat, mark))
}

# ---------------------------------------------------------------------------
# 5. Focus FXN (toujours, même non significatif)
# ---------------------------------------------------------------------------
# FXN peut apparaître dans la colonne 'gene' comme symbole OU comme Ensembl :
# on cherche les deux.
FXN_ENS <- "ENSG00000165060"
fxn_rna <- rna %>% filter(GENE_ID == "FXN" | GENE_ID == FXN_ENS |
                            (!is.na(SYMBOL) & SYMBOL == "FXN") |
                            (!is.na(ENSEMBL) & ENSEMBL == FXN_ENS))

fxn_rows <- list()
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
if (nrow(fxn_rna)) {
  cat(sprintf("RNA-seq : logFC=%.3f  p=%.3g  FDR=%.3g\n",
              fxn_rna$logFC_RNA[1], fxn_rna$PValue_RNA[1], fxn_rna$FDR_RNA[1]))
  write_csv(fxn_rna, file.path(INT_DIR, "FXN_rna.csv"))
} else {
  cat("FXN absent de la table RNA (filtré par min_counts ?).\n")
}

if (nrow(fxn_chip)) {
  print(fxn_chip, row.names = FALSE)
  write_csv(fxn_chip, file.path(INT_DIR, "FXN_summary.csv"))
} else {
  cat("Aucun pic ChIP assigné à FXN (cf. module ciblé 06/07 — c'est attendu\n",
      "pour les marques broad : MACS3 n'appelle pas toujours de pic au locus).\n")
}

cat("\n==> [05b] Terminé. Tables dans results/integration/\n")