#!/usr/bin/env bash
# Decompress IDATs into a work folder (originals are untouched) and verify that
# every sample has BOTH a _Grn and a _Red file.
# usage: bash scripts/01_prepare_idats.sh <raw_dir> <work_dir>
set -euo pipefail
RAW="${1:?raw idat dir}"
WORK="${2:?work dir}"
mkdir -p "$WORK"

n=0
shopt -s nullglob
for f in "$RAW"/*.idat.gz "$RAW"/*.idat; do
  base="$(basename "$f")"; base="${base%.gz}"
  if [ ! -s "$WORK/$base" ]; then
    case "$f" in *.gz) gunzip -c "$f" > "$WORK/$base" ;; *) cp "$f" "$WORK/$base" ;; esac
  fi
  n=$((n + 1))
done
[ "$n" -gt 0 ] || { echo "ERROR: no IDAT files found in the raw folder"; exit 1; }

bad=0
for g in "$WORK"/*_Grn.idat; do
  r="${g%_Grn.idat}_Red.idat"
  [ -s "$r" ] || { echo "ERROR: missing Red file for $(basename "$g")"; bad=1; }
done
for r in "$WORK"/*_Red.idat; do
  g="${r%_Red.idat}_Grn.idat"
  [ -s "$g" ] || { echo "ERROR: missing Grn file for $(basename "$r")"; bad=1; }
done
[ "$bad" -eq 0 ] || exit 1

echo "OK: $(ls "$WORK"/*_Grn.idat | wc -l) arrays (Grn+Red pairs) ready"
