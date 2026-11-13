#!/usr/bin/env bash
##################################################################################################################################
# Step 9: Annotate primary LeGo assembly with Helixer (v0.3.7)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   24_annotation_helixer_gene_models.gff3.gz  Helixer gene models (23,684 genes)
# Input:  $BASE/genomes/lego/lego_pri.fna
# Output: $BASE/helixer/output/lego/lego_pri.helixer.gff3
# Run as a one-task SLURM array on a GPU node, e.g.
#   BASE=$PWD sbatch --array=0 --gres=gpu:1 scripts/09_annotate_genome.sh

set -uo pipefail
die() { echo "ERROR: $*" >&2; exit 1; }

module load miniconda3/3.13
conda activate apptainer
 
BASE=${BASE:?set BASE to the project directory}
SIF=${SIF:-$BASE/helixer/helixer-docker_latest.sif}
GENOMES=($BASE/genomes/lego/lego_pri.fna)
GENOME=${GENOMES[$SLURM_ARRAY_TASK_ID]:-}
 
echo "=== $(date)  node=$(hostname)  task=$SLURM_ARRAY_TASK_ID"
[ -n "$GENOME" ] && [ -s "$GENOME" ] || die "genome not found: '$GENOME'"
[ -e "$SIF" ] || die "container not found: $SIF"
 
NAME=$(basename "$GENOME" .fna)
SP=$(basename "$(dirname "$GENOME")")
OUT=$BASE/helixer/output/$SP
mkdir -p "$OUT" || die "cannot create $OUT"
echo "genome=$GENOME  output=$OUT/$NAME.helixer.gff3"
 
[ -s "$OUT/$NAME.helixer.gff3" ] && { echo "Output exists, skipping."; exit 0; }
 
singularity run --nv --bind "$BASE":"$BASE" --pwd "$BASE" "$SIF" Helixer.py \
  --fasta-path "$GENOME" \
  --lineage vertebrate \
  --gff-output-path "$OUT/$NAME.helixer.gff3.tmp" \
  || die "Helixer exited with an error"
 
N=$(grep -c -P "\tgene\t" "$OUT/$NAME.helixer.gff3.tmp" || true)
[ "${N:-0}" -gt 0 ] || die "output has no gene features"
mv "$OUT/$NAME.helixer.gff3.tmp" "$OUT/$NAME.helixer.gff3"
echo "=== done: $N genes  $(date)"