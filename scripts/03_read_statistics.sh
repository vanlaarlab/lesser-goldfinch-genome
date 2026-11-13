#!/usr/bin/env bash
##################################################################################################################################
# Step 3: Read length and quality statistics of the raw reads (SeqKit v2.1.0)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   03_reads_raw_read_statistics_seqkit.txt  read count, total bases, N50, Q20/Q30 of the raw reads
# Input:  rawlego.fastq (step 2)
# Output: read count, total bases, length quartiles, N50, Q20 and Q30 percentages (printed; Results, Nuclear Genome Assembly)
# Terminal output saved as: seqkit_raw_read_stats.txt
set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
cd "$BASE/assembly"

seqkit stats -a rawlego.fastq | tee seqkit_raw_read_stats.txt
