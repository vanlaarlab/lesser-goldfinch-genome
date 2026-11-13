#!/bin/bash
##################################################################################################################################
# Step 8: Assembly contiguity statistics for both pseudohaplotypes (QUAST v5.2.0)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   05_assembly_contiguity_quast.txt  contiguity statistics, primary and alternate (Table 1)
#   06_assembly_contiguity_quast.tsv  the same as a tab-separated table
# Input:  $BASE/genomes/lego/lego_pri.fna, lego_alt.fna  (GenBank copies of hifiasm hap1 / hap2; step 7)
# Output: $BASE/quast/lego/
# Submit: BASE=$PWD sbatch scripts/08_quast_assemblies.sh

#SBATCH --job-name=quast_lego
#SBATCH --output=logs/quast_lego_%j.out
#SBATCH --error=logs/quast_lego_%j.err
#SBATCH --partition=windfall
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

module load ohpc
module load quast/5.2.0

BASE=${BASE:?set BASE to the project directory}
SP=lego

PRI=$BASE/genomes/$SP/${SP}_pri.fna
ALT=$BASE/genomes/$SP/${SP}_alt.fna

OUT=$BASE/quast/$SP

# Check input files
[ -s "$PRI" ] || { echo "Missing primary assembly: $PRI"; exit 1; }
[ -s "$ALT" ] || { echo "Missing alternate assembly: $ALT"; exit 1; }

mkdir -p "$OUT"

echo "Running QUAST on $SP assemblies..."
echo "Primary: $PRI"
echo "Alternate: $ALT"

quast.py \
    "$PRI" \
    "$ALT" \
    --labels "${SP}_pri,${SP}_alt" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --output-dir "$OUT"

echo "QUAST completed: $(date)"
