#!/usr/bin/env Rscript
# =============================================================================
# 16_integration_diagnostics_plots.R
# -----------------------------------------------------------------------------
# Lit results/integration/integration_*.csv et produit :
#   0. DIAGNOSTIC  — vérifie que les données sont exploitables AVANT de plotter
#   1. Scatter     — logFC ChIP vs logFC RNA, par marque, FXN annoté
#   2. Trend binné — signal ChIP moyen par décile d'expression  <<< le plus fort
#   3. Concordance — test de Fisher sur les signes (robuste au petit n)
#   4. Récap       — table de synthèse
#
# Usage :  Rscript 16_integration_diagnostics_plots.R
#          (depuis la racine du projet, ou export PROJECT_ROOT=/chemin)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel); library(tidyr); library(tibble)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
INT_DIR  <- file.path(PROJECT_ROOT, "results", "integration")
PLOT_DIR <- file.path(INT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

MARK_TYPE <- c(H3K9me3 = "repressive", H3K27me3 = "repressive",
               H3K27ac = "active",     PolII    = "active")

files <- list.files(INT_DIR, pattern = "^integration_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("Aucun integration_*.csv dans ", INT_DIR)

# ---------------------------------------------------------------------------
# Helper : détecte les colonnes quel que soit le suffixe de nom de marque
# ---------------------------------------------------------------------------
# load_one <- function(f) {
#   d <- read_csv(f, show_col_types = FALSE)
#   
#   # nom de la marque = ce qui suit "integration_" et précède le 1er "__" ou ".csv"
#   tag  <- sub("^integration_", "", tools::file_path_sans_ext(basename(f)))
#   mark <- sub("__.*$", "", tag)
#   
#   # colonnes ChIP : logFC_<qqch> / FDR_<qqch> qui ne sont PAS les colonnes RNA
#   lfc_chip_col <- grep("^logFC_", names(d), value = TRUE)
#   lfc_chip_col <- setdiff(lfc_chip_col, "logFC_RNA")
#   fdr_chip_col <- grep("^FDR_", names(d), value = TRUE)
#   fdr_chip_col <- setdiff(fdr_chip_col, "FDR_RNA")
#   
#   if (!length(lfc_chip_col))
#     stop("Pas de colonne logFC ChIP dans ", basename(f),
#          "\nColonnes : ", paste(names(d), collapse = ", "))
#   
#   sym_col <- intersect(c("GeneSymbol", "SYMBOL", "GENE_ID", "gene"), names(d))[1]
#   
#   tibble(
#     mark       = mark,
#     tag        = tag,
#     gene       = d[[sym_col]],
#     logFC_RNA  = suppressWarnings(as.numeric(d$logFC_RNA)),
#     FDR_RNA    = suppressWarnings(as.numeric(d$FDR_RNA)),
#     logFC_chip = suppressWarnings(as.numeric(d[[lfc_chip_col[1]]])),
#     FDR_chip   = suppressWarnings(as.numeric(d[[fdr_chip_col[1]]]))
#   ) %>% filter(!is.na(logFC_RNA), !is.na(logFC_chip))
# }


load_one <- function(f) {
  
  d <- read_csv(f, show_col_types = FALSE)
  
  tag  <- sub("^integration_", "", tools::file_path_sans_ext(basename(f)))
  mark <- sub("__.*$", "", tag)
  
  ## ---------- Gene column ----------
  sym_col <- intersect(
    c("GeneSymbol","GENENAME","gene","SYMBOL","GENE_ID"),
    names(d)
  )[1]
  
  ## ---------- RNA columns ----------
  if ("logFC_RNA" %in% names(d)) {
    logFC_RNA <- d$logFC_RNA
  } else if ("log2FoldChange" %in% names(d)) {
    logFC_RNA <- d$log2FoldChange
  } else {
    stop("No RNA logFC column found in ", basename(f))
  }
  
  if ("FDR_RNA" %in% names(d)) {
    FDR_RNA <- d$FDR_RNA
  } else if ("padj_RNA" %in% names(d)) {
    FDR_RNA <- d$padj_RNA
  } else if ("padj" %in% names(d)) {
    FDR_RNA <- d$padj
  } else {
    FDR_RNA <- rep(NA_real_, nrow(d))
  }
  
  ## ---------- ChIP columns ----------
  lfc_chip_col <- grep("^logFC_", names(d), value = TRUE)
  lfc_chip_col <- setdiff(lfc_chip_col, "logFC_RNA")
  
  fdr_chip_col <- grep("^FDR_", names(d), value = TRUE)
  fdr_chip_col <- setdiff(fdr_chip_col, "FDR_RNA")
  
  if (length(lfc_chip_col) == 0)
    stop("No ChIP logFC column found in ", basename(f))
  
  tibble(
    mark       = mark,
    tag        = tag,
    gene       = d[[sym_col]],
    logFC_RNA  = as.numeric(logFC_RNA),
    FDR_RNA    = as.numeric(FDR_RNA),
    logFC_chip = as.numeric(d[[lfc_chip_col[1]]]),
    FDR_chip   = if (length(fdr_chip_col))
      as.numeric(d[[fdr_chip_col[1]]])
    else
      rep(NA_real_, nrow(d))
  ) %>%
    filter(!is.na(logFC_RNA), !is.na(logFC_chip))
}


dat <- bind_rows(lapply(files, load_one))
dat$type <- MARK_TYPE[dat$mark]

# ===========================================================================
# 0. DIAGNOSTIC — à lire AVANT de regarder les figures
# ===========================================================================
cat("\n=========================== DIAGNOSTIC ===========================\n")
diag <- dat %>%
  group_by(mark) %>%
  summarise(
    n_genes       = n(),
    max_abs_lfc_chip = max(abs(logFC_chip), na.rm = TRUE),
    min_FDR_chip  = min(FDR_chip, na.rm = TRUE),
    sig_chip_05   = sum(FDR_chip < 0.05, na.rm = TRUE),
    max_abs_lfc_rna  = max(abs(logFC_RNA), na.rm = TRUE),
    min_FDR_rna   = min(FDR_RNA, na.rm = TRUE),
    sig_rna_05    = sum(FDR_RNA < 0.05, na.rm = TRUE),
    .groups = "drop")
print(as.data.frame(diag), row.names = FALSE)

cat("\n--- Verdict ---\n")
bad_chip <- diag %>% filter(max_abs_lfc_chip < 0.1)
if (nrow(bad_chip)) {
  cat("!! logFC ChIP ÉCRASÉ pour :", paste(bad_chip$mark, collapse = ", "), "\n")
  cat("   Le fold shrunké de DiffBind 3.x annihile le signal.\n")
  cat("   -> Les scatters seront des bandes verticales. CORRIGE 04 d'abord :\n")
  cat("      logFC <- Conc_FRDA - Conc_Control  (colonnes déjà dans _regions.csv)\n\n")
}
if (all(diag$sig_rna_05 == 0))
  cat("!! Aucun gène RNA significatif (FDR<0.05) : le cross seuillé ne sortira rien.\n",
      "   -> Utilise la CONCORDANCE (section 3), pas les catégories.\n\n")
if (!nrow(bad_chip) && any(diag$sig_rna_05 > 0))
  cat("OK : données exploitables.\n\n")

# ===========================================================================
# 1. SCATTER logFC ChIP vs logFC RNA
# ===========================================================================
# Lecture : marque RÉPRESSIVE -> pente attendue NÉGATIVE (plus de marque = moins d'ARN)
#           marque ACTIVE     -> pente attendue POSITIVE
lab <- dat %>%
  group_by(mark) %>%
  filter(gene == "FXN" | (rank(FDR_RNA, ties.method = "first") <= 8)) %>%
  ungroup()

p1 <- ggplot(dat, aes(logFC_chip, logFC_RNA)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(aes(colour = FDR_RNA < 0.05), alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "#E15759", linewidth = 0.7) +
  geom_text_repel(data = lab, aes(label = gene), size = 2.5,
                  max.overlaps = 20, min.segment.length = 0) +
  geom_point(data = filter(dat, gene == "FXN"),
             colour = "#E15759", size = 3, shape = 21, stroke = 1.2, fill = NA) +
  scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#4E79A7"),
                      name = "RNA FDR < 0.05") +
  facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_x") +
  labs(x = "log2FC ChIP (FRDA vs Control)",
       y = "log2FC RNA (FRDA vs Control)",
       title = "Chromatine vs expression, par marque",
       subtitle = "Répressive : pente attendue négative | Active : pente attendue positive") +
  theme_bw(base_size = 11)

ggsave(file.path(PLOT_DIR, "01_scatter_chip_vs_rna.pdf"), p1, width = 11, height = 8)

# ===========================================================================
# 2. TREND BINNÉ  <<< LE PLUS INFORMATIF À FAIBLE PUISSANCE
# ===========================================================================
# Principe : on range les gènes par déciles de logFC_RNA, et on regarde le
# logFC ChIP MOYEN de chaque décile. Moyenner ~1000 gènes par bin écrase le
# bruit individuel : une tendance monotone apparaît là où le scatter ne montre
# qu'un nuage. C'est la figure standard des papiers ChIP+RNA.
binned <- dat %>%
  group_by(mark) %>%
  mutate(bin = ntile(logFC_RNA, 10)) %>%
  group_by(mark, type, bin) %>%
  summarise(mean_chip = mean(logFC_chip, na.rm = TRUE),
            se_chip   = sd(logFC_chip, na.rm = TRUE) / sqrt(n()),
            n = n(), .groups = "drop")

p2 <- ggplot(binned, aes(bin, mean_chip)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_errorbar(aes(ymin = mean_chip - 1.96 * se_chip,
                    ymax = mean_chip + 1.96 * se_chip),
                width = 0.2, colour = "grey50") +
  geom_line(aes(colour = type), linewidth = 0.8) +
  geom_point(aes(colour = type), size = 2.5) +
  scale_x_continuous(breaks = 1:10) +
  scale_colour_manual(values = c(repressive = "#4E79A7", active = "#E15759")) +
  facet_wrap(~ paste0(mark, " (", type, ")"), scales = "free_y") +
  labs(x = "Décile de log2FC RNA  (1 = plus réprimé  →  10 = plus induit)",
       y = "log2FC ChIP moyen (± IC95)",
       title = "Le signal chromatinien suit-il l'expression ?",
       subtitle = "Répressive : pente descendante attendue | Active : pente montante attendue") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

ggsave(file.path(PLOT_DIR, "02_binned_trend.pdf"), p2, width = 11, height = 8)

# ===========================================================================
# 3. CONCORDANCE DIRECTIONNELLE (Fisher sur les signes)
# ===========================================================================
# Aucun seuil de significativité : on demande juste si les SIGNES s'accordent
# plus souvent que le hasard. Avec des milliers de gènes, ce test a de la
# puissance même quand aucun gène n'est individuellement significatif.
conc <- lapply(split(dat, dat$mark), function(d) {
  d <- d %>% filter(logFC_chip != 0, logFC_RNA != 0)
  tab <- table(chip = sign(d$logFC_chip), rna = sign(d$logFC_RNA))
  if (!all(dim(tab) == c(2, 2)))
    return(tibble(mark = d$mark[1], n = nrow(d), OR = NA, p = NA,
                  pct_concordant = NA, attendu = NA))
  ft <- fisher.test(tab)
  type <- d$type[1]
  # concordant = signes OPPOSÉS si répressif, MÊME signe si actif
  conc_n <- if (type == "repressive") tab["-1", "1"] + tab["1", "-1"]
  else                      tab["-1", "-1"] + tab["1", "1"]
  tibble(mark = d$mark[1], type = type, n = sum(tab),
         OR = round(unname(ft$estimate), 3),
         p  = signif(ft$p.value, 3),
         pct_concordant = round(100 * conc_n / sum(tab), 1),
         attendu = 50)
}) %>% bind_rows()

cat("\n=================== CONCORDANCE DIRECTIONNELLE ===================\n")
cat("(% de gènes dont le signe ChIP s'accorde avec le sens attendu ;\n")
cat(" 50% = hasard. OR et p du test de Fisher sur la table 2x2 des signes.)\n\n")
print(as.data.frame(conc), row.names = FALSE)
write_csv(conc, file.path(INT_DIR, "concordance_summary.csv"))

p3 <- ggplot(conc, aes(reorder(mark, pct_concordant), pct_concordant, fill = type)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey30") +
  geom_text(aes(label = paste0(pct_concordant, "%\np=", p)),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(repressive = "#4E79A7", active = "#E15759")) +
  coord_cartesian(ylim = c(0, 75)) +
  labs(x = NULL, y = "% de gènes concordants avec le sens attendu",
       title = "Concordance chromatine / expression",
       subtitle = "Ligne pointillée = 50% (hasard)") +
  theme_bw(base_size = 11)

ggsave(file.path(PLOT_DIR, "03_concordance.pdf"), p3, width = 7, height = 5)


# ===========================================================================
# 4. BOXPLOTS : distribution du logFC ChIP selon l'état RNA
# ===========================================================================

dat_box <- dat %>%
  mutate(RNA_class = case_when(
    FDR_RNA < 0.05 & logFC_RNA > 0 ~ "Up",
    FDR_RNA < 0.05 & logFC_RNA < 0 ~ "Down",
    TRUE ~ "Stable"
  ))

p4 <- ggplot(dat_box,
             aes(RNA_class, logFC_chip, fill = RNA_class)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~paste0(mark," (",type,")"), scales="free_y") +
  scale_fill_manual(values=c(
    Down="#4E79A7",
    Stable="grey80",
    Up="#E15759")) +
  labs(
    x="RNA category",
    y="ChIP log2FC",
    title="Distribution des changements ChIP selon l'expression RNA") +
  theme_bw(base_size=11) +
  theme(legend.position="none")

ggsave(file.path(PLOT_DIR,
                 "04_chip_boxplot_by_RNAcategory.pdf"),
       p4,width=10,height=7)

library(tidyr)
library(pheatmap)

# ===========================================================================
# 5. HEATMAP des principaux gènes
# ===========================================================================


library(zoo)

roll <- dat %>%
  arrange(logFC_RNA) %>%
  group_by(mark, type) %>%
  mutate(
    rolling = zoo::rollmean(
      logFC_chip,
      k = 250,
      fill = NA,
      align = "center"
    )
  )

p_roll <- ggplot(
  roll,
  aes(logFC_RNA,
      rolling,
      colour=type)
)+
  
  geom_line(linewidth=1)+
  
  facet_wrap(~mark,
             scales="free_y")+
  
  theme_bw()+
  
  labs(
    x="RNA log2FC",
    y="Rolling mean ChIP",
    title="Progressive relationship between chromatin and transcription")


ggsave(
  file.path(PLOT_DIR, "05_rolling_trend.pdf"),
  p_roll,
  width = 11,
  height = 8
)







library(corrplot)



# ===========================================================================
# 7. CORRELATION RNA / ChIP
# ===========================================================================

library(corrplot)

cor_dat <- dat %>%
  group_by(gene, mark) %>%
  summarise(
    logFC_chip = mean(logFC_chip, na.rm=TRUE),
    logFC_RNA  = mean(logFC_RNA, na.rm=TRUE),
    .groups="drop"
  )


# Matrice ChIP
chip_mat <- cor_dat %>%
  dplyr::select(gene, mark, logFC_chip) %>%
  pivot_wider(
    names_from = mark,
    values_from = logFC_chip
  )


# ajouter RNA une seule fois
rna <- cor_dat %>%
  group_by(gene) %>%
  summarise(
    RNA = mean(logFC_RNA, na.rm=TRUE)
  )


chip_mat <- left_join(chip_mat, rna, by="gene")


# transformer en matrice numérique
mat <- chip_mat %>%
  dplyr::select(-gene) %>%
  as.data.frame()

mat <- sapply(mat, as.numeric)

rownames(mat) <- chip_mat$gene


cor_mat <- cor(
  mat,
  use="pairwise.complete.obs",
  method="spearman"
)


pdf(
  file.path(PLOT_DIR,
            "07_correlation_matrix.pdf"),
  width=7,
  height=7
)

corrplot(
  cor_mat,
  method="color",
  addCoef.col="black",
  tl.col="black",
  number.cex=.7
)

dev.off()

# ===========================================================================
# 7. RNA volcano coloré par le ChIP
# ===========================================================================

p5 <- ggplot(dat,
             aes(logFC_RNA,
                 -log10(FDR_RNA),
                 colour=logFC_chip)) +
  
  geom_point(size=.6) +
  
  scale_colour_gradient2(
    low="#4E79A7",
    mid="white",
    high="#E15759") +
  
  facet_wrap(~mark) +
  
  theme_bw(base_size=11) +
  
  labs(
    x="RNA log2FC",
    y="-log10(FDR)",
    colour="ChIP log2FC",
    title="RNA volcano colored by chromatin change")

ggsave(file.path(PLOT_DIR,
                 "07_volcano_chip_colour.pdf"),
       p5,width=11,height=8)

library(pROC)

# ===========================================================================
# 8. ROC : capacité du ChIP à prédire les gènes différentiellement exprimés
# ===========================================================================

roc_dir <- file.path(PLOT_DIR, "ROC")
dir.create(roc_dir, showWarnings = FALSE)

for(m in unique(dat$mark)){
  
  d <- dat %>%
    filter(mark==m)
  
  ## définition des vrais positifs
  d$RNA_sig <- abs(d$logFC_RNA) > 1 & d$FDR_RNA < 0.05
  
  if(length(unique(d$RNA_sig))<2)
    next
  
  rocobj <- roc(
    response=d$RNA_sig,
    predictor=abs(d$logFC_chip),
    quiet=TRUE)
  
  pdf(file.path(roc_dir,
                paste0("ROC_",m,".pdf")),
      width=5,height=5)
  
  plot(
    rocobj,
    legacy.axes=TRUE,
    col="#E15759",
    lwd=3,
    main=paste0(m,
                "\nAUC = ",
                round(auc(rocobj),3)))
  
  abline(a=0,b=1,lty=2,col="grey60")
  
  dev.off()
  
}

auc_table <- lapply(split(dat,dat$mark),function(d){
  
  d$RNA_sig <- abs(d$logFC_RNA)>1 &
    d$FDR_RNA<0.05
  
  if(length(unique(d$RNA_sig))<2)
    return(NULL)
  
  rocobj <- roc(
    d$RNA_sig,
    abs(d$logFC_chip),
    quiet=TRUE)
  
  tibble(
    mark=d$mark[1],
    AUC=as.numeric(auc(rocobj))
  )
  
}) %>%
  bind_rows()

write_csv(
  auc_table,
  file.path(INT_DIR,
            "AUC_summary.csv"))



p_auc <- ggplot(
  auc_table,
  aes(reorder(mark,AUC),AUC,
      fill=AUC>0.7))+
  
  geom_col()+
  
  geom_hline(
    yintercept=.5,
    lty=2)+
  
  coord_flip()+
  
  theme_bw()+
  
  labs(
    x=NULL,
    y="ROC AUC",
    title="Predictive value of chromatin marks")

ggsave(
  file.path(PLOT_DIR,
            "08_AUC_summary.pdf"),
  p_auc,
  width=6,
  height=4)


# ============================================================================
# FXN LOCUS ANALYSIS
# ============================================================================

library(ggplot2)
library(dplyr)

FXN_gene <- "FXN"

fxn_data <- dat %>%
  filter(gene == FXN_gene)


print(fxn_data)


# Scatter FXN dans le contexte global

for(mark_i in unique(dat$mark)){
  
  tmp <- dat %>%
    filter(mark == mark_i)
  
  p <- ggplot(tmp,
              aes(logFC_RNA,
                  logFC_chip))+
    geom_point(alpha=.15)+
    geom_point(
      data=fxn_data %>% filter(mark==mark_i),
      aes(logFC_RNA, logFC_chip),
      size=5
    )+
    geom_hline(yintercept=0,lty=2)+
    geom_vline(xintercept=0,lty=2)+
    ggtitle(
      paste(
        "FXN position -",
        mark_i
      )
    )+
    theme_classic()
  
  
  ggsave(
    file.path(
      PLOT_DIR,
      paste0(
        "FXN_",
        mark_i,
        "_scatter.pdf"
      )
    ),
    p,
    width=6,
    height=5
  )
}


# ============================================================================
# RNA / H3K27ac CONCORDANCE GROUPS
# ============================================================================


h3k27 <- dat %>%
  filter(mark=="H3K27ac") %>%
  dplyr::select(
    gene,
    logFC_RNA,
    logFC_chip
  )


h3k27 <- h3k27 %>%
  mutate(
    group = case_when(
      
      logFC_RNA < 0 &
        logFC_chip < 0
      ~ "RNA_down_H3K27ac_down",
      
      logFC_RNA > 0 &
        logFC_chip > 0
      ~ "RNA_up_H3K27ac_up",
      
      logFC_RNA < 0 &
        logFC_chip > 0
      ~ "Discordant_down",
      
      logFC_RNA > 0 &
        logFC_chip < 0
      ~ "Discordant_up"
      
    )
  )


table(h3k27$group)



write.csv(
  h3k27,
  file.path(
    PLOT_DIR,
    "H3K27ac_RNA_concordance_groups.csv"
  ),
  row.names=FALSE
)



install.packages("clusterProfiler")
install.packages("org.Hs.eg.db")


# ============================================================================
# GO ENRICHMENT
# ============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)


genes_down <- h3k27 %>%
  filter(
    group=="RNA_down_H3K27ac_down"
  ) %>%
  pull(gene)



ego <- enrichGO(
  gene = genes_down,
  OrgDb = org.Hs.eg.db,
  keyType="SYMBOL",
  ont="BP",
  pAdjustMethod="BH"
)


pdf(
  file.path(
    PLOT_DIR,
    "GO_RNA_down_H3K27ac_down.pdf"
  ),
  width=8,
  height=6
)


dotplot(
  ego,
  showCategory=20
)


dev.off()




install.packages("msigdbr")


# ============================================================================
# MITOCHONDRIAL ENRICHMENT
# ============================================================================


library(msigdbr)
library(dplyr)


# Récupération des gènes GO mitochondriaux

mito_genes <- msigdbr(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
) %>%
  filter(
    grepl(
      "mitochond",
      gs_name,
      ignore.case = TRUE
    )
  ) %>%
  pull(gene_symbol)


length(mito_genes)


# ============================================================================
# INTEGRATIVE HEATMAP
# ============================================================================


library(pheatmap)


heat <- dat %>%
  filter(
    mark %in%
      c(
        "H3K27ac",
        "PolII"
      )
  ) %>%
  dplyr::select(
    gene,
    mark,
    logFC_chip
  ) %>%
  pivot_wider(
    names_from=mark,
    values_from=logFC_chip
  )


heat <- left_join(
  heat,
  dat %>%
    dplyr::select(
      gene,
      logFC_RNA
    ) %>%
    distinct(),
  by="gene"
)



top_genes <- heat %>%
  arrange(
    abs(logFC_RNA),
    decreasing=TRUE
  ) %>%
  head(100)



mat <- top_genes %>%
  column_to_rownames("gene") %>%
  as.matrix()



pdf(
  file.path(
    PLOT_DIR,
    "Heatmap_RNA_H3K27ac_PolII.pdf"
  ),
  width=7,
  height=10
)

pheatmap(
  mat,
  scale="row"
)

dev.off()




# ============================================================================
# REGION TYPE ANALYSIS
# ============================================================================


dat %>%
  group_by(
    mark,
    region
  )%>%
  summarise(
    mean_change=
      mean(
        logFC_chip,
        na.rm=TRUE
      )
  )

ggplot(
  region_summary,
  aes(
    region,
    mean_change,
    fill=mark
  )
)+
  geom_col(position="dodge")+
  theme_classic()





ggplot(
  dat,
  aes(
    distance_TSS,
    logFC_chip,
    colour=mark
  )
)+
  geom_smooth()+
  theme_classic()











# ===========================================================================
# 4. RÉCAP + FXN
# ===========================================================================
cat("\n=========================== FXN ===========================\n")
fxn <- dat %>% filter(gene == "FXN") %>%
  dplyr::select(mark, type, logFC_chip, FDR_chip, logFC_RNA, FDR_RNA)
if (nrow(fxn)) {
  print(as.data.frame(fxn), row.names = FALSE)
  write_csv(fxn, file.path(INT_DIR, "FXN_cross.csv"))
} else {
  cat("FXN absent des tables d'intégration.\n")
  cat("-> Normal pour les marques broad : aucun pic MACS3 assigné au locus.\n")
  cat("-> Utilise le module ciblé (06/07 + fxn_summary_CvsF.csv) pour FXN.\n")
}

cat("\n==> Figures : ", PLOT_DIR, "\n")
cat("    01_scatter_chip_vs_rna.pdf\n")
cat("    02_binned_trend.pdf        <-- commence par celle-ci\n")
cat("    03_concordance.pdf\n\n")