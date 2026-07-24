#!/usr/bin/env bash
# =============================================================================
# 08_fxn_bigwig.sh — Pistes de couverture (bigWig) pour IGV, focus FXN
# -----------------------------------------------------------------------------
# Génère des bigWig normalisés (CPM) restreints au locus FXN -> légers et
# instantanés à ouvrir dans IGV (génome hg38).
# Réindexe d'abord les BAM pour supprimer l'avertissement "index older".
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
require_cmd bamCoverage
require_cmd samtools

OUTDIR="$RESULTS_DIR/fxn/bigwig"
mkdir -p "$OUTDIR"
REGION="chr9:68900000:69200000"   # FXN + marges

echo "==> [08] Réindexation des BAM (corrige l'avertissement d'index)"
for b in "$BAM_DIR"/*.bam; do samtools index "$b" 2>/dev/null || true; done

echo "==> [08] Génération des bigWig (région $REGION)"
for b in "$BAM_DIR"/*.bam; do
  name="$(basename "$b" .bam)"
  bamCoverage -b "$b" -o "$OUTDIR/${name}.bw" \
    --normalizeUsing CPM --binSize 10 --region "$REGION" -p "$THREADS" \
    2>/dev/null
  echo "    $name.bw"
done

cat <<MSG

==> [08] Terminé. bigWig dans : $OUTDIR

Pour visualiser :
  1. Ouvre IGV :  igv     (sinon : sudo apt install igv)
  2. Genome -> Human (hg38)
  3. File -> Load from File -> sélectionne les .bw
  4. Tape 'FXN' dans la barre de recherche.

À regarder :
  - H3K9me3 : signal plus haut chez les F (FRDA) sur l'intron 1 / corps du gène ?
  - H3K27ac & PolII : signal plus bas chez les F au promoteur ?
  - Y a-t-il un TROU de couverture sur la zone (GAA)n (intron 1) dans TOUS les
    échantillons, input compris ? -> ce serait un artefact de mappabilité.
MSG
