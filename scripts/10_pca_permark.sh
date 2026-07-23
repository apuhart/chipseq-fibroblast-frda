#!/usr/bin/env bash
# =============================================================================
# 10_pca_permark.sh — PCA & corrélation par marque (diagnostic des "0")
# -----------------------------------------------------------------------------
# Pour chaque marque, calcule un signal par échantillon sur les régions
# consensus (>=2 éch.), puis produit une PCA et une heatmap de corrélation.
# But : voir si les échantillons se séparent par CONDITION (Control vs FRDA)
# ou par autre chose (batch, variabilité individuelle) -> explique pourquoi
# peu de régions passent le FDR.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
require_cmd multiBamSummary
require_cmd plotPCA
require_cmd plotCorrelation

OUTDIR="$RESULTS_DIR/figures/pca"
mkdir -p "$OUTDIR"

for mark in H3K27ac H3K27me3 H3K9me3 PolII; do
  bed="$CONSENSUS_DIR/${mark}_consensus_min2.bed"
  [[ -s "$bed" ]] || { echo "[WARN] $bed absent, marque ignorée."; continue; }

  # BAM + labels de cette marque (depuis le manifest)
  mapfile -t bams < <(awk -F, -v m="$mark" '$2==m{print $4}' "$MANIFEST_PEAKS")
  labels=$(for b in "${bams[@]}"; do basename "$b" | cut -d- -f1; done | tr '\n' ' ')
  echo "==> [$mark] ${#bams[@]} échantillons"

  npz="$OUTDIR/${mark}.npz"
  multiBamSummary BED-file --BED "$bed" \
    --bamfiles "${bams[@]}" --labels $labels \
    -p "$THREADS" -o "$npz" 2> "$LOG_DIR/pca_${mark}.log"

  plotPCA -in "$npz" -o "$OUTDIR/${mark}_PCA.pdf" \
    -T "PCA — $mark" 2>> "$LOG_DIR/pca_${mark}.log"

  plotCorrelation -in "$npz" --corMethod spearman --whatToPlot heatmap \
    --plotNumbers -o "$OUTDIR/${mark}_corr.pdf" \
    -T "Corrélation — $mark" 2>> "$LOG_DIR/pca_${mark}.log"

  echo "    -> ${mark}_PCA.pdf  ${mark}_corr.pdf"
done

cat <<'MSG'

==> [10] Terminé. Figures dans results/figures/pca/

Lecture :
  - Si sur la PCA les Control et les FRDA se séparent nettement le long de
    PC1/PC2 -> le signal différentiel existe (les "0" viennent alors surtout
    du seuil FDR strict avec n=3).
  - S'ils se mélangent -> la variabilité individuelle domine l'effet maladie :
    c'est la vraie raison du faible nombre de régions significatives.
MSG
