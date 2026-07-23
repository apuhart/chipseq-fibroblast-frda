#!/usr/bin/env bash
# =============================================================================
# run_all.sh — Orchestrateur des étapes bash (00 -> 03)
# -----------------------------------------------------------------------------
# Les étapes R (04, 05) sont lancées séparément car DiffBind nécessite que tu
# aies complété metadata/samples_metadata.csv (conditions, réplicats).
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash "$HERE/scripts/00_setup.sh"

cat <<'MSG'

------------------------------------------------------------------
ACTION REQUISE : édite metadata/samples_metadata.csv
  - colonne Condition : ex. FRDA / Control / Edited
  - colonne Replicate : entier unique par groupe (>=2 par condition)
Puis appuie sur Entrée pour continuer (Ctrl-C pour arrêter).
------------------------------------------------------------------
MSG
read -r _

# bash "$HERE/scripts/01_qc.sh"
# bash "$HERE/scripts/02_macs3_peakcalling.sh"
# bash "$HERE/scripts/03_consensus.sh"
# Rscript "$HERE/scripts/04_diffbind_chipseeker.R"
# Rscript "$HERE/scripts/05b_integrate_edgeR.R"
# Rscript "$HERE/scripts/05_integrate_rnaseq.R"
# Rscript "$HERE/scripts/05b_integrate_deseq2.R"
# bash "$HERE/scripts/03_consensus.sh"
# Rscript "$HERE/scripts/05b_integrate_edgeR.R"
# bash "$HERE/scripts/06b_fxn_targeted_custom.sh"
# Rscript "$HERE/scripts/07_fxn_analysis.R"
# bash "$HERE/scripts/08_fxn_bigwig.sh"
# Rscript "$HERE/scripts/09_volcano.R"
# bash "$HERE/scripts/10_pca_permark.sh"
# Rscript "$HERE/scripts/11_integration_plots.R"
# Rscript "$HERE/scripts/12_fxn_synthesis.R"
# Rscript "$HERE/scripts/13_pathway_enrichment.R"
# Rscript "$HERE/scripts/14_multimark_concordance.R"
# Rscript "$HERE/scripts/15_fxn_neighborhood.R"



Rscript "$HERE/scripts/04b_recompute_gene_level.R"
# Rscript "$HERE/scripts/05b_integrate_deseq2.R"
Rscript "$HERE/scripts/16_integration_diagnostics_plots3.R"








cat <<'MSG'

==> Pipeline terminées.

MSG
