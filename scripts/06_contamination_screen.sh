#!/usr/bin/env bash
##################################################################################################################################
# Step 6: Screen for adapter and foreign DNA contamination, contig level (NCBI FCS v0.5.5: FCS-adaptor and FCS-GX)
##################################################################################################################################
# Supporting files: none
# Input:  legoContigs_hap1.fasta, legoContigs_hap2.fasta (step 5)
# Output: outputdirHap1/, outputdirHap2/ (FCS-adaptor); gx_outLeGoHap1/, gx_outLeGoHap2/ (FCS-GX)
# run_fcsadaptor.sh, fcs.py and the FCS-GX database (db/) are expected in $BASE/assembly; see
# https://github.com/ncbi/fcs for downloads. Tax ID 54559 = Spinus psaltria (NCBI Taxonomy).
# On an HPC without Docker, add --container-engine singularity --image fcs-adaptor.sif to
# run_fcsadaptor.sh and --image fcs-gx.sif to fcs.py.
set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
mkdir -p "$BASE/assembly"
cd "$BASE/assembly"

# Adaptor contamination screen
./run_fcsadaptor.sh --fasta-input legoContigs_hap1.fasta --output-dir ./outputdirHap1 --euk
./run_fcsadaptor.sh --fasta-input legoContigs_hap2.fasta --output-dir ./outputdirHap2 --euk

# Foreign DNA contamination screen
python3 ./fcs.py screen genome --fasta legoContigs_hap1.fasta --out-dir ./gx_outLeGoHap1/ --gx-db "$PWD/db" --tax-id 54559
python3 ./fcs.py screen genome --fasta legoContigs_hap2.fasta --out-dir ./gx_outLeGoHap2/ --gx-db "$PWD/db" --tax-id 54559
