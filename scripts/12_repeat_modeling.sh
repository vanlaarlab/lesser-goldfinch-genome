#!/usr/bin/env bash
##################################################################################################################################
# Step 12: Identify repeat families de novo on the primary assembly (RepeatModeler v2.0.6)
##################################################################################################################################
# Supporting files: none (the raw families are curated in step 13)
# Input:  legoContigs_hap1.fasta (hap1 / primary contigs from step 5)
# Output: $BASE/repeats/legoRepeats-families.fa (raw RepeatModeler2 library; curated in step 13)
# Runs for several days; keep it in a screen/tmux session or submit it as a job.
set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
mkdir -p "$BASE/repeats"
cd "$BASE/repeats"
ln -sf "$BASE/assembly/legoContigs_hap1.fasta" .

BuildDatabase -name legoRepeats legoContigs_hap1.fasta
RepeatModeler -database legoRepeats -threads 12 -LTRStruct >& run.out
