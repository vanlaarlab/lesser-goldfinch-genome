#!/usr/bin/env bash
##################################################################################################################################
# Step 4: Remove reads that still have adapters (cutadapt v5.0)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   02_reads_adapter_trimming_cutadapt.txt  Cutadapt report (2 reads with adapters removed)
# Input:  rawlego.fastq
# Output: trimmedlego.fastq  (reads with a detected adapter, in either
#         orientation, are discarded rather than trimmed)
# Terminal output saved as: 2_removeAdapters.txt
set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
mkdir -p "$BASE/assembly"
cd "$BASE/assembly"

cutadapt \
  -a First_adapter=ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT \
  -a Second_adapter=ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT \
  --error-rate=0.1 \
  --overlap=35 \
  --revcomp \
  --discard-trimmed \
  -o trimmedlego.fastq \
  rawlego.fastq
