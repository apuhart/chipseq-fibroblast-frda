#!/usr/bin/env bash
# =============================================================================
# run_fibroblasts_all.sh — Pipeline fibroblastes complet, de bout en bout
# -----------------------------------------------------------------------------
# Enchaîne toutes les étapes dans l'ordre. À lancer depuis la RACINE du projet,
# environnement conda/venv activé (macs3, samtools, bedtools, deeptools, R).
#
#   cd /home/apuhart/Downloads/pipeline_chiprna
#   source ~/chipseq-env/bin/activate
#   bash scripts/run_fibroblasts_all.sh
#
# Prérequis : metadata/samples_metadata.csv DÉJÀ rempli (colonne Condition).
# Le script s'arrête à la moindre erreur et journalise tout dans logs/.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_ROOT="$(cd "$HERE/.." && pwd)"
source "$HERE/config.sh"
mkdir -p "$LOG_DIR"

# Journal global horodaté
RUN_LOG="$LOG_DIR/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$RUN_LOG") 2>&1
echo "================ Pipeline fibroblastes — $(date) ================"
echo "PROJECT_ROOT = $PROJECT_ROOT"
echo "Journal      = $RUN_LOG"

# Petite fonction d'étape : annonce + exécute + vérifie
step() {
  local title="$1"; shift
  echo ""
  echo "---------------------------------------------------------------"
  echo ">>> $title"
  echo "---------------------------------------------------------------"
  "$@"
}

# --- 00 : détection des échantillons (n'écrase pas un metadata existant) ---
step "[00] Setup / détection des échantillons" bash "$HERE/00_setup.sh"

# --- Garde-fou : la colonne Condition doit être remplie ---
if awk -F, 'NR>1 && ($3=="" || $3=="NA")' "$METADATA" | grep -q .; then
  echo ""
  echo "[ARRÊT] La colonne 'Condition' de $METADATA n'est pas entièrement remplie."
  echo "        Complète-la (Control / FRDA ...) puis relance ce script."
  exit 1
fi
echo "    Métadonnées OK :"
cut -d, -f1,3,5 "$METADATA" | sed 's/^/      /'

# --- Étapes principales ---
step "[01] Contrôle qualité"                 bash    "$HERE/01_qc.sh"
step "[02] Appel de pics MACS3"              bash    "$HERE/02_macs3_peakcalling.sh"
step "[03] Régions consensus"                bash    "$HERE/03_consensus.sh"
step "[04] DiffBind + ChIPseeker"            Rscript "$HERE/04_diffbind_chipseeker.R"

# --- Intégration RNA-seq (edgeR / Ensembl) ---
step "[05b] Intégration RNA-seq (edgeR)"     Rscript "$HERE/05b_integrate_edgeR.R"

# --- Analyse ciblée FXN ---
step "[06] Comptage ciblé FXN"               bash    "$HERE/06_fxn_targeted.sh"
step "[07] Analyse quantitative FXN"         Rscript "$HERE/07_fxn_analysis.R"
step "[08] bigWig FXN pour IGV"              bash    "$HERE/08_fxn_bigwig.sh"

# --- Diagnostics visuels ---
step "[09] Volcano par marque"               Rscript "$HERE/09_volcano.R"
step "[10] PCA & corrélation par marque"     bash    "$HERE/10_pca_permark.sh"

# --- Récapitulatif ---
echo ""
echo "================ TERMINÉ — $(date) ================"
echo "Résultats principaux :"
echo "  QC ............... results/qc/qc_summary.csv"
echo "  Pics ............. results/peaks/<MARK>/"
echo "  Différentiel ..... results/diffbind/*_regions.csv"
echo "  Niveau gène ...... results/annotation/*_gene_level.csv"
echo "  Intégration RNA .. results/integration/integration_*.csv"
echo "  FXN (ciblé) ...... results/fxn/fxn_summary_CvsF.csv , fxn_barplots.pdf"
echo "  FXN (bigWig IGV) . results/fxn/bigwig/*.bw"
echo "  Figures .......... results/figures/ (volcano, pca/)"
echo "  Journal complet .. $RUN_LOG"
