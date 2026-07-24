#!/usr/bin/env bash
# =============================================================================
# 02_macs3_peakcalling.sh — Appel de pics avec MACS3
# -----------------------------------------------------------------------------
# Lit metadata/chip_manifest.csv et, pour chaque ChIP :
#   - détecte automatiquement paired-end vs single-end (format BAMPE / BAM) ;
#   - utilise l'input apparié comme contrôle (-c) s'il existe ;
#   - lance MACS3 en mode broad (H3K9me3, H3K27me3) ou narrow (H3K27ac, PolII) ;
#   - écrit un manifest enrichi (chemins des fichiers de pics) pour l'étape R.
# Sorties : results/peaks/<MARK>/<SAMPLE>_<MARK>_peaks.{broad|narrow}Peak
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
require_cmd macs3
require_cmd samtools

[ -f "$MANIFEST" ] || { echo "[ERREUR] Manifest absent. Lance d'abord 00_setup.sh." >&2; exit 1; }

echo "==> [02] Appel de pics MACS3"
echo "sample,mark,peaktype,chip_bam,input_bam,peaks" > "$MANIFEST_PEAKS"

# Détection PE/SE robuste (sans déclencher SIGPIPE sous 'set -o pipefail')
is_paired() {
  local bam="$1" first
  set +o pipefail
  first="$(samtools view -f 1 "$bam" 2>/dev/null | head -1)"
  set -o pipefail
  [[ -n "$first" ]]
}

# On saute l'en-tête du manifest et on boucle sur chaque ligne.
tail -n +2 "$MANIFEST" | while IFS=, read -r sample mark ptype chip input; do
  [[ -z "$sample" ]] && continue
  outdir="$PEAKS_DIR/$mark"
  mkdir -p "$outdir"
  name="${sample}_${mark}"

  # Format d'entrée
  if is_paired "$chip"; then fmt="BAMPE"; else fmt="BAM"; fi

  # Argument de contrôle (input) seulement s'il existe
  ctrl_args=()
  if [[ -n "$input" && -f "$input" ]]; then
    ctrl_args=(-c "$input")
  else
    echo "    [WARN] $name : appel SANS input (résultats moins fiables)."
  fi

  echo "    $name  [type=$ptype | format=$fmt]"

  if [[ "$ptype" == "broad" ]]; then
    macs3 callpeak \
      -t "$chip" "${ctrl_args[@]}" \
      -f "$fmt" -g "$GENOME_SIZE" \
      --broad --broad-cutoff "$BROAD_CUTOFF" \
      -q "$QVALUE" --keep-dup "$KEEPDUP" \
      -n "$name" --outdir "$outdir" \
      &> "$LOG_DIR/macs3_${name}.log"
    peaks="$outdir/${name}_peaks.broadPeak"
  else
    macs3 callpeak \
      -t "$chip" "${ctrl_args[@]}" \
      -f "$fmt" -g "$GENOME_SIZE" \
      -q "$QVALUE" --keep-dup "$KEEPDUP" \
      -n "$name" --outdir "$outdir" \
      &> "$LOG_DIR/macs3_${name}.log"
    peaks="$outdir/${name}_peaks.narrowPeak"
  fi

  if [[ ! -s "$peaks" ]]; then
    echo "    [WARN] $name : aucun pic produit (voir $LOG_DIR/macs3_${name}.log)."
  fi
  echo "$sample,$mark,$ptype,$chip,$input,$peaks" >> "$MANIFEST_PEAKS"
done

echo "    Manifest enrichi : $MANIFEST_PEAKS"
echo "==> [02] Terminé."
