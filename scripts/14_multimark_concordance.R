#!/usr/bin/env Rscript
# =============================================================================
# 14_multimark_concordance.R — Gènes concordants sur plusieurs marques
# -----------------------------------------------------------------------------
# Croise les 4 tables d'intégration + l'expression RNA pour repérer les gènes
# où PLUSIEURS marques bougent de façon cohérente avec un silencing :
#   gain H3K9me3 (+), gain H3K27me3 (+), perte H3K27ac (-), perte PolII (-),
#   baisse d'expression RNA (-).
# Un gène soutenu par plusieurs marques est bien plus crédible qu'un hit isolé
# -> filtre de robustesse face au manque de puissance (n=3).
#
# Sortie : results/integration/multimark_concordance.csv (trié par score)
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(purrr) })

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
INT <- file.path(PROJECT_ROOT, "results", "integration")

read_mark <- function(mark) {
  f <- file.path(INT, paste0("integration_", mark, ".csv"))
  if (!file.exists(f)) return(NULL)
  d <- read_csv(f, show_col_types = FALSE)
  lfc <- paste0("logFC_", mark)
  d %>% transmute(ENTREZID = as.character(ENTREZID), GeneSymbol,
                  logFC_RNA, FDR_RNA,
                  !!paste0("lfc_", mark) := .data[[lfc]])
}

marks <- c("H3K9me3", "H3K27me3", "H3K27ac", "PolII")
tabs <- compact(map(marks, read_mark))
if (length(tabs) == 0) stop("Aucune table d'intégration trouvée.")

# Jointure progressive sur ENTREZID (garde RNA depuis la 1re table)
base <- tabs[[1]]
for (i in seq_along(tabs)[-1]) {
  add <- tabs[[i]] %>% select(ENTREZID, starts_with("lfc_"))
  base <- full_join(base, add, by = "ENTREZID")
}

# Score de concordance "silencing" : +1 par marque allant dans le sens attendu
score_dir <- function(x, expect_sign) ifelse(is.na(x), 0, as.integer(sign(x) == expect_sign))
res <- base %>%
  mutate(
    s_H3K9me3  = if ("lfc_H3K9me3"  %in% names(.)) score_dir(lfc_H3K9me3,  +1) else 0,
    s_H3K27me3 = if ("lfc_H3K27me3" %in% names(.)) score_dir(lfc_H3K27me3, +1) else 0,
    s_H3K27ac  = if ("lfc_H3K27ac"  %in% names(.)) score_dir(lfc_H3K27ac,  -1) else 0,
    s_PolII    = if ("lfc_PolII"    %in% names(.)) score_dir(lfc_PolII,    -1) else 0,
    s_RNA      = score_dir(logFC_RNA, -1),
    concordance_score = s_H3K9me3 + s_H3K27me3 + s_H3K27ac + s_PolII + s_RNA
  ) %>%
  arrange(desc(concordance_score), logFC_RNA)

write_csv(res, file.path(INT, "multimark_concordance.csv"))
cat("Concordance -> ", file.path(INT, "multimark_concordance.csv"), "\n\n")

cat("Distribution du score de concordance (0 = aucun sens attendu, 5 = tout concorde) :\n")
print(table(res$concordance_score))
cat("\nTop 20 gènes les plus concordants avec un silencing épigénétique :\n")
print(as.data.frame(head(res %>%
  select(GeneSymbol, logFC_RNA, starts_with("lfc_"), concordance_score), 20)),
  row.names = FALSE)
cat("\nNOTE : exploratoire. Un score élevé = hypothèse à valider, pas une preuve.\n")
