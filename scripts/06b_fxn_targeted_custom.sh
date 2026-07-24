#!/usr/bin/env bash
# =============================================================================
# 06b_fxn_targeted_custom.sh — Comptage ciblé FXN sur un BED EXISTANT
# -----------------------------------------------------------------------------
# Identique à 06_fxn_targeted.sh MAIS n'écrase PAS le fichier de régions :
# il utilise tel quel results/fxn/fxn_regions.bed (que tu as édité à la main).
# Permet de tester différents découpages (flancs de l'expansion GAA, etc.)
# sans retoucher au reste du pipeline.
#
# Sortie : results/fxn/fxn_counts.csv  (puis Rscript scripts/07_fxn_analysis.R)
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
require_cmd samtools

OUTDIR="$RESULTS_DIR/fxn"
BED="$OUTDIR/fxn_regions.bed"

[[ -s "$BED" ]] || { echo "[ERREUR] $BED introuvable ou vide. Crée-le d'abord." >&2; exit 1; }

echo "==> [06b] Comptage ciblé FXN (BED fourni, NON écrasé)"
echo "    Régions utilisées :"; sed 's/^/      /' "$BED"

declare -A LIBSIZE
libsize() {
  local bam="$1"
  if [[ -z "${LIBSIZE[$bam]:-}" ]]; then
    LIBSIZE[$bam]=$(samtools idxstats "$bam" | awk '{s+=$3} END{print s}')
  fi
  echo "${LIBSIZE[$bam]}"
}
count_region() { samtools view -c -F 260 "$1" "$2" 2>/dev/null; }

OUT="$OUTDIR/fxn_counts.csv"
echo "sample,mark,region,chip_count,chip_libsize,input_count,input_libsize" > "$OUT"

tail -n +2 "$MANIFEST_PEAKS" | while IFS=, read -r sample mark ptype chip input peaks; do
  [[ -z "$sample" ]] && continue
  clib=$(libsize "$chip")
  ilib=""; [[ -n "$input" && -f "$input" ]] && ilib=$(libsize "$input")
  while read -r chr start end name; do
    [[ -z "$chr" || "$chr" == \#* ]] && continue
    region="${chr}:${start}-${end}"
    cc=$(count_region "$chip" "$region")
    ic=""; [[ -n "$input" && -f "$input" ]] && ic=$(count_region "$input" "$region")
    echo "$sample,$mark,$name,$cc,$clib,${ic:-NA},${ilib:-NA}" >> "$OUT"
  done < "$BED"
done

echo "    Table -> $OUT"
echo "    Mappabilité (input) par sous-région :"
awk -F, 'NR>1 {key=$3; sum[key]+=$6; n[key]++}
         END{for(k in sum) printf "      %-22s input moyen = %.0f reads\n", k, sum[k]/n[k]}' "$OUT"

echo "==> [06b] Terminé. Lance : Rscript scripts/07_fxn_analysis.R"
