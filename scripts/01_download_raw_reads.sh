#!/usr/bin/env bash
##################################################################################################################################
# Step 1: Download raw PacBio HiFi reads from UC Davis (BioShare)
##################################################################################################################################
# Supporting files: none (raw reads are in NCBI SRA, BioSample SAMN49551328)
# Credentials are NOT stored in this repository. Before running, export:
#   BIOSHARE_URL         full wget link to hifi_reads/wget_index.html sent by UC Davis
#   BIOSHARE_SESSIONID   your BioShare session id (browser cookie)
# Raw reads are also archived at NCBI (BioSample SAMN49551328; see README).
#
# Output: *.bam HiFi reads (barcoded) in $BASE/assembly
set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
mkdir -p "$BASE/assembly"
cd "$BASE/assembly"
: "${BIOSHARE_URL:?set BIOSHARE_URL}"
: "${BIOSHARE_SESSIONID:?set BIOSHARE_SESSIONID}"

wget -r --level=10 -nH -nc --cut-dirs=3 --no-parent \
  --reject "wget_index.html" --no-check-certificate \
  --header "Cookie: sessionid=${BIOSHARE_SESSIONID};" \
  "$BIOSHARE_URL"
