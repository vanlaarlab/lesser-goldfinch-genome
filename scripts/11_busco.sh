#!/bin/bash
##################################################################################################################################
# Step 11: BUSCO completeness of both assemblies and the predicted proteins (BUSCO v6.1.0)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   07_assembly_busco_primary.txt     BUSCO, primary assembly (97.8%)
#   08_assembly_busco_alternate.txt   BUSCO, alternate assembly (95.5%)
#   27_annotation_busco_proteins.txt  BUSCO, predicted proteins (94.1%)
# =============================================================================
# BUSCO for the lesser goldfinch with the passeriformes odb12.2 dataset
#
#   task 0  lego_proteins  protein mode  eggnog/proteins/lego_pri.faa (Helixer; from step 10)
#   task 1  lego_pri       genome mode   genomes/lego/lego_pri.fna
#   task 2  lego_alt       genome mode   genomes/lego/lego_alt.fna
#
# Usage, from the project directory (export BASE=$PWD):
#   1. once, on the login node (needs internet; downloads the dataset):
#        bash scripts/11_busco.sh setup
#   2. submit the three runs (they work offline from the downloaded dataset);
#      set cpus, memory and time for your cluster:
#        sbatch --array=0-2 --cpus-per-task=<N> --mem=<M> --time=<T> scripts/11_busco.sh
#   3. when all three have finished:
#        bash scripts/11_busco.sh summary
#
# Requeue-safe: a finished task is skipped; an interrupted one restarts with
# BUSCO --restart. Outputs go to busco/<name>/ and busco/busco_summary.tsv.
# Environment: conda env "busco" (BUSCO 6.1.0).
# =============================================================================

set -euo pipefail

BASE=${BASE:?set BASE to the project directory}
OUT=$BASE/busco
DL=$OUT/busco_downloads
LINEAGE_FILE=$OUT/lineage.txt          # the exact dataset name chosen at setup

NAMES=(lego_proteins lego_pri lego_alt)
INPUTS=("$BASE/eggnog/proteins/lego_pri.faa"
        "$BASE/genomes/lego/lego_pri.fna"
        "$BASE/genomes/lego/lego_alt.fna")
MODES=(protein genome genome)

die() { echo "ERROR: $*" >&2; exit 1; }

module load miniconda3/3.13
eval "$(conda shell.bash hook)"
conda activate busco

# -----------------------------------------------------------------------------
# setup: choose and download the dataset once
# -----------------------------------------------------------------------------
if [ "${1:-}" = "setup" ]; then
  mkdir -p "$DL" "$BASE/logs"
  busco --version
  for i in 0 1 2; do
    [ -s "${INPUTS[$i]}" ] || die "missing input: ${INPUTS[$i]}"
  done

  avail=$(busco --list-datasets 2>&1 | grep -o 'passeriformes_odb1[0-9][.0-9]*' | sort -u || true)
  echo "passeriformes datasets offered: ${avail:-none}"
  if   grep -qx 'passeriformes_odb12.2' <<< "$avail"; then L=passeriformes_odb12.2
  elif grep -qx 'passeriformes_odb12'   <<< "$avail"; then L=passeriformes_odb12
       echo "WARNING: odb12.2 not offered; using odb12 (a different dataset)" >&2
  else die "no passeriformes odb12 dataset offered by this BUSCO version"
  fi

  busco --download "$L" --download_path "$DL"
  [ -d "$DL/lineages/$L" ] || die "download did not create $DL/lineages/$L"
  echo "$L" > "$LINEAGE_FILE"
  echo "dataset: $L ($(ls "$DL/lineages/$L/hmms" | wc -l) BUSCO profiles)"
  echo "next: sbatch --array=0-2 ... scripts/11_busco.sh"
  exit 0
fi

# -----------------------------------------------------------------------------
# summary: one row per run
# -----------------------------------------------------------------------------
if [ "${1:-}" = "summary" ]; then
  L=$(cat "$LINEAGE_FILE" 2>/dev/null) || die "run setup first"
  {
    printf "run\tmode\tdataset\tcomplete\tsingle\tduplicated\tfragmented\tmissing\tn\n"
    for i in 0 1 2; do
      f=$(ls "$OUT/${NAMES[$i]}"/short_summary*.txt 2>/dev/null | head -1 || true)
      if [ -z "$f" ]; then
        printf "%s\t%s\t%s\tnot finished\n" "${NAMES[$i]}" "${MODES[$i]}" "$L"
        continue
      fi
      grep -m1 -oE 'C:[0-9.]+%\[S:[0-9.]+%,D:[0-9.]+%\],F:[0-9.]+%,M:[0-9.]+%,n:[0-9]+' "$f" \
        | sed -E 's/C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%,n:([0-9]+)/\1\t\2\t\3\t\4\t\5\t\6/' \
        | sed "s/^/${NAMES[$i]}\t${MODES[$i]}\t$L\t/"
    done
  } | tee "$OUT/busco_summary.tsv"
  exit 0
fi

# -----------------------------------------------------------------------------
# run: one array task per input
# -----------------------------------------------------------------------------
i=${SLURM_ARRAY_TASK_ID:?"submit with sbatch, or use 'setup' / 'summary'"}
NAME=${NAMES[$i]}; IN=${INPUTS[$i]}; MODE=${MODES[$i]}
L=$(cat "$LINEAGE_FILE" 2>/dev/null) || die "run 'bash scripts/11_busco.sh setup' first"
[ -d "$DL/lineages/$L" ] || die "dataset not downloaded: $DL/lineages/$L"
[ -s "$IN" ] || die "input not found on $(hostname): $IN"

echo "task $i on $(hostname): $NAME ($MODE) with $L"
busco --version

if ls "$OUT/$NAME"/short_summary*.txt >/dev/null 2>&1; then
  echo "already finished, skipping"
  exit 0
fi
RESTART=""
[ -d "$OUT/$NAME" ] && RESTART="--restart" && echo "resuming previous run"

cd "$OUT"
busco -i "$IN" -m "$MODE" \
  -l "$DL/lineages/$L" \
  -o "$NAME" --out_path "$OUT" \
  -c "$SLURM_CPUS_PER_TASK" \
  --offline --download_path "$DL" \
  $RESTART

grep -m1 -E 'C:[0-9.]+%' "$OUT/$NAME"/short_summary*.txt
