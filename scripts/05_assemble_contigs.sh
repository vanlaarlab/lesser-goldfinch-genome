#!/usr/bin/env bash
##################################################################################################################################
# Step 5: Assemble nuclear genome, contig level (hifiasm v0.25.0)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   04_assembly_hifiasm_log.txt  hifiasm log: version, command, coverage peaks
# Input:  trimmedlego.fastq
# Output: legoContigs_hap1.fasta (hap1 / primary pseudohaplotype)
#         legoContigs_hap2.fasta (hap2 / alternate pseudohaplotype)
# GFA -> FASTA conversion is done with awk (no scaffolding step is used).
# Terminal output saved as: 4_assembleNuclearContigs.txt
set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
mkdir -p "$BASE/assembly"
cd "$BASE/assembly"

hifiasm -o legoContigs -t 16 trimmedlego.fastq
awk '/^S/{print ">"$2"\n"$3}' legoContigs.bp.hap1.p_ctg.gfa | fold > legoContigs_hap1.fasta
awk '/^S/{print ">"$2"\n"$3}' legoContigs.bp.hap2.p_ctg.gfa | fold > legoContigs_hap2.fasta
