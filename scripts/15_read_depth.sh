#!/usr/bin/env bash
##################################################################################################################################
# Step 15: Read depth per scaffold, to check the sex-chromosome assignment (NCBI datasets v18.28.0; minimap2 v2.31; SAMtools v1.24)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   12_sexchr_read_depth_samtools_coverage.tsv  per-scaffold read depth (Z/W depth check)
# Input:  $BASE/assembly/trimmedlego.fastq (step 4); sexchr/lego_pri.sexchr_assignment.tsv (step 14)
# Output: $BASE/depth/lego_coverage.tsv (per-scaffold depth) and the length-weighted mean depth of autosomal, Z and W
#         scaffolds (Results: 34.8x, 19.1x and 18.0x; in a female, Z and W are expected at about half the autosomal depth)
set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
mkdir -p "$BASE/depth"
cd "$BASE/depth"
ln -sf "$BASE/assembly/trimmedlego.fastq" trimmedlego.fastq

# Step 15a: Download LEGO hap1 reference (NCBI datasets v18.28.0)

datasets download genome accession GCA_060154605.1 --include genome,gff3
unzip ncbi_dataset.zip
mv ncbi_dataset/data/GCA_060154605.1/*.fna lesser_goldfinch.fa

# Step 15b: Build split index for low-memory mapping (minimap2 v2.31)

minimap2 -x map-hifi -I 500M -d lesser_goldfinch.mmi lesser_goldfinch.fa

# Step 15c: Map trimmed HiFi reads to hap1 (minimap2 v2.31, samtools v1.24)

minimap2 -ax map-hifi -t 3 -K 50M --split-prefix lego_tmp \
  lesser_goldfinch.mmi \
  trimmedlego.fastq \
  | samtools view -b -o lego.unsorted.bam -

# Step 15d: Sort and index alignments (samtools v1.24)

samtools quickcheck lego.unsorted.bam && echo "OK"
samtools sort -@ 3 -m 768M -o lego.sorted.bam lego.unsorted.bam
samtools index -@ 3 lego.sorted.bam

# Step 15e: Mapping statistics (samtools v1.24)

samtools flagstat lego.sorted.bam

##################################################################
# Step 15f: Per-scaffold read depth (samtools v1.24)

samtools coverage lego.sorted.bam > lego_coverage.tsv

##################################################################
# Step 15g: Length-weighted mean depth by sex-chromosome class (class = column 11 of the step 14 assignment table)

awk -F'\t' 'NR == FNR { if (FNR > 1) c[$1] = $11; next }
     FNR > 1 { k = ($1 in c) ? c[$1] : "NA"; L = $3 - $2 + 1; len[k] += L; d[k] += L * $7; n[k]++ }
     END { for (k in len) printf "%-10s %4d scaffolds %9.1f Mb  length-weighted mean depth %5.1fx\n", k, n[k], len[k] / 1e6, d[k] / len[k] }' \
  "$BASE/sexchr/lego_pri.sexchr_assignment.tsv" lego_coverage.tsv | sort
