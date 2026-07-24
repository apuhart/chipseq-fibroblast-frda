#!/usr/bin/env bash
# =============================================================================
# 01_qc.sh — Contrôle qualité basique des BAM
# -----------------------------------------------------------------------------
# - Indexe les BAM si nécessaire (.bai requis par idxstats et deeptools/DiffBind).
# - samtools flagstat  : taux d'alignement, duplicats, reads appariés.
# - samtools idxstats  : reads par chromosome.
# - (optionnel) deeptools : plotFingerprint + corrélation inter-échantillons,
#   exécuté seulement si deeptools est installé.
# Toutes les sorties vont dans results/qc/.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
require_cmd samtools

echo "==> [01] Contrôle qualité"
SUMMARY="$QC_DIR/qc_summary.csv"
echo "bam,total_reads,mapped_reads,mapped_pct,duplicates,properly_paired" > "$SUMMARY"

index_if_needed() {
  local bam="$1"
  if [ ! -f "${bam}.bai" ] && [ ! -f "${bam%.bam}.bai" ]; then
    # Essaie d'indexer à côté du BAM ; si le dossier est en lecture seule,
    # écrit l'index dans results/qc/.
    if ! samtools index -@ "$THREADS" "$bam" 2>/dev/null; then
      local alt="$QC_DIR/$(basename "$bam").bai"
      echo "    [INFO] Index écrit dans $alt (dossier source non inscriptible)."
      samtools index -@ "$THREADS" "$bam" "$alt"
    fi
  fi
}

shopt -s nullglob
for bam in "$BAM_DIR"/*"$BAM_SUFFIX"; do
  name="$(basename "$bam" "$BAM_SUFFIX")"
  echo "    QC : $name"

  index_if_needed "$bam" || echo "    [WARN] Indexation impossible pour $name."

  fs="$QC_DIR/${name}.flagstat.txt"
  samtools flagstat -@ "$THREADS" "$bam" > "$fs"
  samtools idxstats "$bam" > "$QC_DIR/${name}.idxstats.txt" 2>/dev/null || true

  # Extraction de quelques chiffres clés vers le récapitulatif CSV
  total=$(grep -m1 "in total"        "$fs" | awk '{print $1}')
  mapped=$(grep -m1 " mapped ("      "$fs" | awk '{print $1}')
  mapped_pct=$(grep -m1 " mapped ("  "$fs" | sed -E 's/.*\(([0-9.]+%).*/\1/')
  dups=$(grep -m1 "duplicates"       "$fs" | awk '{print $1}')
  paired=$(grep -m1 "properly paired" "$fs" | awk '{print $1}')
  echo "$name,${total:-NA},${mapped:-NA},${mapped_pct:-NA},${dups:-NA},${paired:-NA}" >> "$SUMMARY"
done
shopt -u nullglob
echo "    Récapitulatif : $SUMMARY"

# ---------------------------------------------------------------------------
# QC avancé optionnel (deeptools)
# ---------------------------------------------------------------------------
if command -v multiBamSummary >/dev/null 2>&1 && command -v plotCorrelation >/dev/null 2>&1; then
  echo "    deeptools détecté : corrélation inter-échantillons + fingerprint."
  mapfile -t ALLBAMS < <(ls "$BAM_DIR"/*"$BAM_SUFFIX")
  LABELS=$(for b in "${ALLBAMS[@]}"; do basename "$b" "$BAM_SUFFIX"; done | tr '\n' ' ')

  multiBamSummary bins \
    --bamfiles "${ALLBAMS[@]}" \
    --labels $LABELS \
    --numberOfProcessors "$THREADS" \
    -o "$QC_DIR/multiBamSummary.npz" 2> "$LOG_DIR/multiBamSummary.log" || \
    echo "    [WARN] multiBamSummary a échoué (voir log)."

  if [ -f "$QC_DIR/multiBamSummary.npz" ]; then
    plotCorrelation -in "$QC_DIR/multiBamSummary.npz" \
      --corMethod spearman --whatToPlot heatmap \
      --plotFile "$QC_DIR/sample_correlation.pdf" \
      --outFileCorMatrix "$QC_DIR/sample_correlation_matrix.tsv" \
      2>> "$LOG_DIR/multiBamSummary.log" || true
  fi

  if command -v plotFingerprint >/dev/null 2>&1; then
    plotFingerprint --bamfiles "${ALLBAMS[@]}" --labels $LABELS \
      --numberOfProcessors "$THREADS" \
      --plotFile "$QC_DIR/fingerprint.pdf" \
      2> "$LOG_DIR/fingerprint.log" || true
  fi
else
  echo "    (deeptools non installé : QC avancé ignoré — facultatif)"
fi

echo "==> [01] Terminé."
