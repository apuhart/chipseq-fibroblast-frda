#!/usr/bin/env bash
# =============================================================================
# 00_setup.sh
# -----------------------------------------------------------------------------
# 1. Crée l'arborescence de résultats.
# 2. Détecte automatiquement les échantillons à partir des noms de BAM.
# 3. Associe chaque ChIP à l'input correspondant (même préfixe SAMPLE).
# 4. Écrit un manifest (metadata/chip_manifest.csv).
# 5. Génère un TEMPLATE de métadonnées à compléter par l'utilisateur
#    (conditions, traitements, réplicats) -> indispensable pour DiffBind.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"

echo "==> [00] Création de l'arborescence"
mkdir -p "$QC_DIR" "$PEAKS_DIR" "$CONSENSUS_DIR" "$LOG_DIR" "$META_DIR" \
         "$RESULTS_DIR/diffbind" "$RESULTS_DIR/annotation" "$RESULTS_DIR/integration"

# Vérifie qu'il y a bien des BAM
shopt -s nullglob
BAMS=("$BAM_DIR"/*"$BAM_SUFFIX")
shopt -u nullglob
if [ "${#BAMS[@]}" -eq 0 ]; then
  echo "[ERREUR] Aucun fichier *$BAM_SUFFIX dans $BAM_DIR" >&2
  exit 1
fi
echo "    ${#BAMS[@]} fichiers BAM détectés."

# Fonction d'extraction des champs depuis un nom de fichier
parse_sample() { local f; f="$(basename "$1")"; echo "${f%%-*}"; }                 # avant le 1er '-'
parse_mark()   { local f r; f="$(basename "$1")"; r="${f#*-}"; echo "${r%%_*}"; }  # entre '-' et '_'

# ---------------------------------------------------------------------------
# Passe 1 : recenser les inputs par échantillon
# ---------------------------------------------------------------------------
declare -A INPUT_BY_SAMPLE
for bam in "${BAMS[@]}"; do
  sample="$(parse_sample "$bam")"
  mark="$(parse_mark "$bam")"
  # comparaison insensible à la casse pour le tag input
  if [[ "${mark,,}" == "${INPUT_TAG,,}" ]]; then
    INPUT_BY_SAMPLE["$sample"]="$bam"
  fi
done

# ---------------------------------------------------------------------------
# Passe 2 : écrire le manifest ChIP (en excluant les inputs)
# ---------------------------------------------------------------------------
echo "sample,mark,peaktype,chip_bam,input_bam" > "$MANIFEST"
declare -A SEEN_SAMPLES
for bam in "${BAMS[@]}"; do
  sample="$(parse_sample "$bam")"
  mark="$(parse_mark "$bam")"
  [[ "${mark,,}" == "${INPUT_TAG,,}" ]] && continue

  # type de pic (défaut narrow si marque non listée)
  ptype="${PEAKTYPE[$mark]:-narrow}"
  if [[ -z "${PEAKTYPE[$mark]:-}" ]]; then
    echo "    [WARN] Marque inconnue '$mark' -> traitée en narrow par défaut."
  fi

  inbam="${INPUT_BY_SAMPLE[$sample]:-}"
  if [[ -z "$inbam" ]]; then
    # Échantillon incomplet (pas d'input) -> on l'ÉCARTE du pipeline.
    echo "    [SKIP] '$sample' ($mark) : aucun input -> échantillon ignoré."
    continue
  fi

  SEEN_SAMPLES["$sample"]=1   # uniquement les échantillons retenus (avec input)
  echo "$sample,$mark,$ptype,$bam,$inbam" >> "$MANIFEST"
done
echo "    Manifest écrit : $MANIFEST"

# ---------------------------------------------------------------------------
# Template de métadonnées : 1 ligne par couple (sample x mark) pour permettre
# des réplicats biologiques propres. À COMPLÉTER avant l'étape DiffBind.
# ---------------------------------------------------------------------------
if [ ! -f "$METADATA" ]; then
  echo "sample,Patient,Condition,Treatment,Replicate" > "$METADATA"
  for s in "${!SEEN_SAMPLES[@]}"; do
    # Valeurs par défaut : à éditer (Condition vide = à remplir).
    echo "$s,$s,,unedited,1" >> "$METADATA"
  done
  echo "    Template de métadonnées créé : $METADATA"
  echo "    >>> ÉDITE-LE : remplis 'Condition' (ex: FRDA / Control / Edited)"
  echo "    >>> et 'Replicate' (entier unique par groupe) avant DiffBind."
else
  echo "    Métadonnées déjà présentes (non écrasées) : $METADATA"
fi

echo "==> [00] Terminé."
