#!/usr/bin/env bash
##################################################################################################################################
# Step 2: Convert BAM to FASTQ (samtools v1.21)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   01_reads_bam_to_fastq_samtools.txt  BAM to FASTQ conversion log (2,939,006 reads)
# Input:  m84066_250331_024525_s4.hifi_reads.bc2089.bam
# Output: rawlego.fastq
# Terminal output saved as: 1_convertBAMtoFASTQ.txt
set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
mkdir -p "$BASE/assembly"
cd "$BASE/assembly"

samtools fastq m84066_250331_024525_s4.hifi_reads.bc2089.bam > rawlego.fastq
