#!/usr/bin/env bash
# =============================================================================
# 06_fxn_targeted.sh — Comptage ciblé sur le locus FXN (hg38)
# -----------------------------------------------------------------------------
# Contourne l'analyse par pics (aveugle à FXN) en comptant DIRECTEMENT les
# lectures sur des sous-régions de FXN, pour chaque marque, ChIP ET input.
# Sortie : results/fxn/fxn_counts.csv  (exploité par 07_fxn_analysis.R)
#
# Coordonnées FXN hg38 (GRCh38), brin + :
#   gène      chr9:69,035,563-69,100,178
#   la région (GAA)n pathogène est dans l'INTRON 1 (~chr9:69,037,000-69,041,000)
# Ajuste les bornes ci-dessous si besoin.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
require_cmd samtools

OUTDIR="$RESULTS_DIR/fxn"
mkdir -p "$OUTDIR"
BED="$OUTDIR/fxn_regions.bed"

# Sous-régions d'intérêt (chrom  start  end  nom)
cat > "$BED" <<'EOF'
chr9	69034563	69036563	FXN_promoter
chr9	69037000	69041000	FXN_intron1_GAA
chr9	69035563	69100178	FXN_genebody
EOF
echo "==> [06] Comptage ciblé FXN"
echo "    Sous-régions :"; sed 's/^/      /' "$BED"

# Taille de librairie (lectures mappées) via l'index : instantané.
declare -A LIBSIZE
libsize() {
  local bam="$1"
  if [[ -z "${LIBSIZE[$bam]:-}" ]]; then
    LIBSIZE[$bam]=$(samtools idxstats "$bam" | awk '{s+=$3} END{print s}')
  fi
  echo "${LIBSIZE[$bam]}"
}

# Comptage de lectures mappées primaires (-F 260) dans une région.
count_region() { samtools view -c -F 260 "$1" "$2" 2>/dev/null; }

OUT="$OUTDIR/fxn_counts.csv"
echo "sample,mark,region,chip_count,chip_libsize,input_count,input_libsize" > "$OUT"

tail -n +2 "$MANIFEST_PEAKS" | while IFS=, read -r sample mark ptype chip input peaks; do
  [[ -z "$sample" ]] && continue
  clib=$(libsize "$chip")
  ilib=""; [[ -n "$input" && -f "$input" ]] && ilib=$(libsize "$input")

  while read -r chr start end name; do
    region="${chr}:${start}-${end}"
    cc=$(count_region "$chip" "$region")
    ic=""; [[ -n "$input" && -f "$input" ]] && ic=$(count_region "$input" "$region")
    echo "$sample,$mark,$name,$cc,$clib,${ic:-NA},${ilib:-NA}" >> "$OUT"
  done < "$BED"
done

echo "    Table de comptage -> $OUT"

# --- Détection d'un éventuel trou de couverture sur la région (GAA)n ---
# Si même l'input a ~0 lecture sur intron1_GAA, c'est un problème de
# mappabilité (séquence répétée) : FXN serait alors invisible pour raison
# technique, pas biologique.
echo "    Vérification mappabilité région GAA (input) :"
awk -F, 'NR>1 && $3=="FXN_intron1_GAA" {printf "      %-8s input=%s reads\n", $1, $6}' "$OUT"

echo "==> [06] Terminé. Lance ensuite : Rscript scripts/07_fxn_analysis.R"
