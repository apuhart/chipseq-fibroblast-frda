# Pipeline ChIP-seq — projet Ataxie de Friedreich (fibroblastes / iNeurons)

Préparation des données ChIP-seq et intégration avec le RNA-seq (DEGs).

## Arborescence

```
frda_chipseq_pipeline/
├── scripts/
│   ├── config.sh                  # chemins, génome, seuils (À ÉDITER)
│   ├── 00_setup.sh                # dossiers + détection échantillons + manifest + template métadonnées
│   ├── 01_qc.sh                   # QC samtools (+ deeptools si dispo)
│   ├── 02_macs3_peakcalling.sh    # MACS3 broad/narrow, PE/SE auto, input apparié
│   ├── 03_consensus.sh            # régions consensus par marque (bedtools)
│   ├── 04_diffbind_chipseeker.R   # DiffBind + ChIPseeker + matrice niveau gène
│   ├── 05_integrate_rnaseq.R      # fusion ChIP <-> RNA-seq + classification
│   └── run_all.sh                 # orchestrateur étapes bash
├── metadata/                      # généré : manifests + samples_metadata.csv (À COMPLÉTER)
├── results/
│   ├── qc/                        # flagstat, idxstats, corrélations
│   ├── peaks/<MARK>/              # *.narrowPeak / *.broadPeak
│   ├── consensus/                 # *_consensus_union.bed, *_consensus_min2.bed
│   ├── diffbind/                  # <MARK>__<CONTRASTE>_regions.csv
│   ├── annotation/                # <MARK>__<CONTRASTE>_gene_level.csv, *_consensus_annotation.csv
│   └── integration/               # integration_*.csv, candidates_*.csv
├── logs/                          # logs MACS3 et deeptools
└── rnaseq/DEGs.csv                # À FOURNIR : ta table de DEGs RNA-seq
```

## Dépendances

- **bash** : `macs3`, `samtools`, `bedtools` (et `deeptools`, optionnel).
- **R / Bioconductor** : `DiffBind`, `ChIPseeker`, `TxDb.Hsapiens.UCSC.hg38.knownGene`,
  `org.Hs.eg.db`, `GenomicRanges`, `GenomeInfoDb`, plus `dplyr`, `readr`, `tidyr`, `stringr`, `purrr`.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DiffBind","ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db"))
install.packages(c("dplyr","readr","tidyr","stringr","purrr"))
```

## Utilisation

```bash
cd frda_chipseq_pipeline
# 1) Vérifie/édite scripts/config.sh (BAM_DIR, génome, seuils)
bash scripts/run_all.sh          # 00 -> 03 (s'arrête pour les métadonnées)

# 2) Édite metadata/samples_metadata.csv : Condition + Replicate (>=2 par groupe)

# 3) Étapes R
export PROJECT_ROOT="$(pwd)"
Rscript scripts/04_diffbind_chipseeker.R
# adapte RNA_FILE / colonnes dans scripts/05_integrate_rnaseq.R puis :
Rscript scripts/05_integrate_rnaseq.R
```

## Table finale

`results/integration/integration_<MARK>__<CONTRASTE>.csv` contient :

| GeneID | GeneSymbol | logFC_RNA | padj_RNA | signal_<MARK> | logFC_<MARK> | FDR_<MARK> | category |
|--------|------------|-----------|----------|---------------|--------------|------------|----------|

`category` repère les profils attendus (ex. `gain_mark_down_expr` pour H3K9me3
gagné + expression réduite). `candidates_*.csv` ne garde que ces profils.

## Points d'attention

- **Génome hg38** supposé : si autre build, change `GENOME_BUILD` (config.sh) et le `TxDb` (04).
- **Réplicats** : DiffBind exige ≥2 réplicats par condition et par marque, sinon
  l'étape différentielle est sautée (annotation seule).
- **Style de chromosomes** harmonisé en UCSC (`chr*`) côté R pour matcher le TxDb.
- **PolII en narrow** et **H3K9me3 via MACS3** : choix respectés mais discutables
  (PolII souvent broad ; H3K9me3 parfois mieux capté par epic2/SICER).
- **Agrégation pic→gène** : par défaut, région la plus significative par gène
  (`AGGREGATION <- "best_fdr"` dans 04, modifiable en `"max_abs"`).
- **Fusion sur GeneSymbol** : les `GeneID` (Entrez) sont conservés pour lever les
  ambiguïtés de symboles si besoin.
