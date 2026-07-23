#!/usr/bin/env bash
# =============================================================================
# config.sh — Configuration centrale du pipeline ChIP-seq (projet FRDA)
# -----------------------------------------------------------------------------
# Ce fichier est "sourcé" par tous les autres scripts bash.
# Modifie ici les chemins, le génome et les seuils : rien d'autre à toucher.
# =============================================================================

# Racine du projet = dossier parent de scripts/
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ---------------------------------------------------------------------------
# ENTREES
# ---------------------------------------------------------------------------
# Dossier contenant les BAM déjà alignés et triés.
export BAM_DIR="/home/apuhart/Downloads/ChIPseq/DATA/BAM/Fibroblasts"

# Convention de nommage attendue : {SAMPLE}-{MARK}_{TAG}_sorted.bam
#   ex : C6718-H3K9me3_MN2_sorted.bam  -> SAMPLE=C6718, MARK=H3K9me3
#        C6718-input_MN1_sorted.bam    -> input du SAMPLE C6718
export BAM_SUFFIX="_sorted.bam"   # suffixe commun à tous les BAM
export INPUT_TAG="input"          # mot-clé identifiant les inputs (insensible casse plus bas)

# ---------------------------------------------------------------------------
# GENOME
# ---------------------------------------------------------------------------
# Taille effective du génome pour MACS3 : hs (humain), mm (souris) ou un nombre.
export GENOME_SIZE="hs"
# Build (purement informatif côté bash ; côté R il faut le TxDb correspondant).
export GENOME_BUILD="hg38"

# ---------------------------------------------------------------------------
# SEUILS MACS3
# ---------------------------------------------------------------------------
export QVALUE="0.05"          # seuil q-value
export BROAD_CUTOFF="0.1"     # seuil pour les régions broad
export KEEPDUP="auto"         # gestion des duplicats (auto | all | 1)

# ---------------------------------------------------------------------------
# MARQUES -> TYPE DE PIC
# (tableau associatif : disponible car config.sh est sourcé dans le même shell)
# ---------------------------------------------------------------------------
declare -gA PEAKTYPE=(
  [H3K9me3]=broad
  [H3K27me3]=broad
  [H3K27ac]=narrow
  [PolII]=narrow
)
# Marque inconnue -> traitée en narrow par défaut (voir 00_setup.sh).

# ---------------------------------------------------------------------------
# SORTIES
# ---------------------------------------------------------------------------
export RESULTS_DIR="$PROJECT_ROOT/results"
export PEAKS_DIR="$RESULTS_DIR/peaks"
export QC_DIR="$RESULTS_DIR/qc"
export CONSENSUS_DIR="$RESULTS_DIR/consensus"
export LOG_DIR="$PROJECT_ROOT/logs"
export META_DIR="$PROJECT_ROOT/metadata"

export THREADS="${THREADS:-4}"

# Fichiers manifest (générés automatiquement)
export MANIFEST="$META_DIR/chip_manifest.csv"
export MANIFEST_PEAKS="$META_DIR/chip_manifest_with_peaks.csv"
export METADATA="$META_DIR/samples_metadata.csv"

# ---------------------------------------------------------------------------
# Petit utilitaire : vérifie qu'un binaire est installé.
# ---------------------------------------------------------------------------
require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "[ERREUR] '$1' introuvable dans le PATH. Installe-le avant de continuer." >&2
    return 1
  fi
}
