#!/usr/bin/env bash
# =============================================================================
# 03_consensus.sh — Régions consensus par marque
# -----------------------------------------------------------------------------
# Pour chaque marque, construit :
#   - <MARK>_consensus_union.bed : union de tous les pics (merge).
#   - <MARK>_consensus_min2.bed  : régions présentes dans >=2 échantillons
#                                  (seulement s'il y a >=2 fichiers de pics).
# Remarque : DiffBind reconstruit de toute façon son propre consensus à partir
# des BAM ; ces BED servent de référence descriptive (étape 4 du cahier des
# charges) et pour des visualisations.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
require_cmd bedtools

[ -f "$MANIFEST_PEAKS" ] || { echo "[ERREUR] Lance d'abord 02_macs3_peakcalling.sh." >&2; exit 1; }

echo "==> [03] Construction des consensus par marque"

# Liste des marques présentes dans le manifest enrichi
marks=$(tail -n +2 "$MANIFEST_PEAKS" | cut -d, -f2 | sort -u)

for mark in $marks; do
  # Fichiers de pics existants pour cette marque
  mapfile -t files < <(tail -n +2 "$MANIFEST_PEAKS" \
                        | awk -F, -v m="$mark" '$2==m && $6!="" {print $6}' \
                        | while read -r p; do [ -s "$p" ] && echo "$p"; done)

  n="${#files[@]}"
  if [ "$n" -eq 0 ]; then
    echo "    [WARN] $mark : aucun fichier de pics non vide, consensus ignoré."
    continue
  fi
  echo "    $mark : $n fichier(s) de pics"

  union="$CONSENSUS_DIR/${mark}_consensus_union.bed"
  cat "${files[@]}" | cut -f1-3 | sort -k1,1 -k2,2n \
    | bedtools merge -i - > "$union"
  echo "      union  -> $union ($(wc -l < "$union") régions)"

  if [ "$n" -ge 2 ]; then
    min2="$CONSENSUS_DIR/${mark}_consensus_min2.bed"
    # multiinter trouve les segments couverts par >=1 fichier ; colonne 4 = nb fichiers
    bedtools multiinter -i "${files[@]}" \
      | awk 'BEGIN{OFS="\t"} $4>=2 {print $1,$2,$3}' \
      | sort -k1,1 -k2,2n | bedtools merge -i - > "$min2"
    echo "      >=2 éch -> $min2 ($(wc -l < "$min2") régions)"
  fi
done

echo "==> [03] Terminé."
