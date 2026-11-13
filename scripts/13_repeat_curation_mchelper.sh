#!/usr/bin/env bash
##################################################################################################################################
# Step 13: Curate the repeat library (MCHelper v1.7.1, Docker image built from the MCHelper repository)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   13_repeats_curated_library.fa                 curated repeat library (368 families; deposited on Dryad)
#   15_repeats_library_curation_mchelper_log.txt  MCHelper run log
# Input:  $BASE/repeats/legoRepeats-families.fa (step 12), legoContigs_hap1.fasta (step 5),
#         BUSCO passeriformes_odb12 HMM profiles
# Output: $BASE/repeats/mchelp_out/curated_sequences_NR.fa (curated library, 368 families; deposited on Dryad)
#         copied to $BASE/te_lego/library/lego_curated_sequences_NR.fa for step 14
# MCHelper settings: automatic mode (-a F), classified + unclassified modules (-r A), 1,000 bp extension (-e 1000).
set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
cd "$BASE/repeats"

# 1. MCHelper container (once)
[ -d MCHelper ] || git clone https://github.com/GonzalezLab/MCHelper.git
docker build -t mchelper MCHelper

# 2. BUSCO HMM profiles used to remove gene-derived families: the passeriformes_odb12 dataset's profiles
#    concatenated into one file (this preparation step was not captured in the original terminal log)
busco --download passeriformes_odb12 --download_path busco_downloads
cat busco_downloads/lineages/passeriformes_odb12/hmms/*.hmm > passeriformes_odb12.hmm

# 3. Pfam profiles for MCHelper's database, then curation
docker run --rm -v "$PWD":/work mchelper /bin/bash -c '
  cd /opt/mchelper/db &&
  wget -q https://urgi.versailles.inrae.fr/download/repet/profiles/ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz &&
  tar xf ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz &&
  mv ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm Pfam35.0.hmm &&
  cd /work &&
  python3 /opt/mchelper/MCHelper.py -l legoRepeats-families.fa -g legoContigs_hap1.fasta -o mchelp_out \
    -t 16 -e 1000 --input_type fasta -a F -r A -b passeriformes_odb12.hmm'

mkdir -p "$BASE/te_lego/library"
cp mchelp_out/curated_sequences_NR.fa "$BASE/te_lego/library/lego_curated_sequences_NR.fa"
grep -c '>' "$BASE/te_lego/library/lego_curated_sequences_NR.fa"
