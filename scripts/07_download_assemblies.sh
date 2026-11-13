#!/usr/bin/env bash
##################################################################################################################################
# Step 7: Download the assemblies and reference annotations used by steps 8-19 (NCBI datasets v18.28.0)
##################################################################################################################################
# Supporting files: none (public NCBI data)
# The lesser goldfinch assemblies analyzed in steps 8-19 are the GenBank copies of the hifiasm contigs from step 5
# (identical sequence; GenBank contig names JBZEFH...). Also downloaded: the other finches (TLR analysis; the house
# finch also serves as the reference for sex-chromosome assignment), zebra finch, chicken and human.
#
# Output layout (relative to $BASE):
#   genomes/lego/lego_pri.fna (+ .fai), genomes/lego/lego_alt.fna         lesser goldfinch, primary and alternate
#   genomes/hofi/hofi_pri.fna, genomes/pisi/pisi_pri.fna, genomes/amgo/amgo_pri.fna
#   genomes/tagu/bTaeGut7_mat.fna and genomes/tagu/ncbi_dataset/data/GCF_048771995.1/{protein.faa,cds_from_genomic.fna,genomic.gff}
#   gaga_prot.faa, gaga_genomic.gff, refs/chicken/cds_from_genomic.fna     chicken
#   refs/human/genomic.gff.gz                                              human T2T-CHM13
set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
cd "$BASE"
TMP=$(mktemp -d)

# get <accession> <files to include> -> unpacked into $TMP/<accession>
get() {
  datasets download genome accession "$1" --include "$2" --filename "$TMP/$1.zip"
  unzip -q -o "$TMP/$1.zip" -d "$TMP/$1"
}
genomic() { find "$TMP/$1" -name "*_genomic.fna" | head -1; }

mkdir -p genomes/lego genomes/hofi genomes/pisi genomes/amgo genomes/tagu refs/chicken refs/human

# lesser goldfinch: primary (GCA_060154605.1) and alternate (BioProject PRJNA1281196)
get GCA_060154605.1 genome && cp "$(genomic GCA_060154605.1)" genomes/lego/lego_pri.fna
get PRJNA1281196    genome && cp "$(genomic PRJNA1281196)"    genomes/lego/lego_alt.fna
samtools faidx genomes/lego/lego_pri.fna

# other finches (primary assemblies; GenBank headers label the Z and W chromosomes)
get GCA_027477595.1 genome && cp "$(genomic GCA_027477595.1)" genomes/hofi/hofi_pri.fna     # house finch, bHaeMex1
get GCA_053477405.1 genome && cp "$(genomic GCA_053477405.1)" genomes/pisi/pisi_pri.fna     # pine siskin
get GCA_057381415.1 genome && cp "$(genomic GCA_057381415.1)" genomes/amgo/amgo_pri.fna     # American goldfinch, bSpiTri1

# zebra finch: GenBank genome (bTaeGut7 maternal) and RefSeq annotation
get GCA_048771995.1 genome && cp "$(genomic GCA_048771995.1)" genomes/tagu/bTaeGut7_mat.fna
get GCF_048771995.1 protein,cds,gff3
mkdir -p genomes/tagu/ncbi_dataset/data/GCF_048771995.1
cp "$TMP"/GCF_048771995.1/ncbi_dataset/data/GCF_048771995.1/{protein.faa,cds_from_genomic.fna,genomic.gff} \
   genomes/tagu/ncbi_dataset/data/GCF_048771995.1/

# chicken (GRCg7w, GCF_016700215.2): proteins, annotation and CDS
get GCF_016700215.2 protein,cds,gff3
D=$TMP/GCF_016700215.2/ncbi_dataset/data/GCF_016700215.2
cp "$D/protein.faa" gaga_prot.faa
cp "$D/genomic.gff" gaga_genomic.gff
cp "$D/cds_from_genomic.fna" refs/chicken/cds_from_genomic.fna

# human T2T-CHM13v2.0 annotation (TLR1-family comparison)
get GCF_009914755.1 gff3
gzip -c "$TMP/GCF_009914755.1/ncbi_dataset/data/GCF_009914755.1/genomic.gff" > refs/human/genomic.gff.gz

rm -rf "$TMP"
for f in genomes/*/*.fna; do echo "$f: $(grep -c '>' "$f") sequences"; done
