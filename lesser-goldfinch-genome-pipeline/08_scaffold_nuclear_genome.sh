# Step 7: Scaffolding nuclear genome (scaffold-level) (ntLink v1.3.11 / ragtag)

gunzip bSpiTri1.HiC.hap1.20240426.fasta.gz
gunzip bSpiTri1.HiC.hap2.20240426.fasta.gz

# BLAST to identify CHD1W sex chromosome locus
makeblastdb -in bSpiTri1.HiC.hap1.20240426.fasta -dbtype nucl -out hap1_db
makeblastdb -in bSpiTri1.HiC.hap2.20240426.fasta -dbtype nucl -out hap2_db

blastn -query CHD1W.fa -db hap1_db -out CHD1W_vs_hap1.tsv \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

blastn -query CHD1W.fa -db hap2_db -out CHD1W_vs_hap2.tsv \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# Scaffold hap1 contigs onto reference haplotypes
ragtag.py scaffold -u bSpiTri1.HiC.hap1.20240426.fasta legoContigs_hap1.fasta -o ragtag_legoHap1_on_AmGoHap1
ragtag.py scaffold -u bSpiTri1.HiC.hap2.20240426.fasta legoContigs_hap1.fasta -o ragtag_legoHap1_on_AmGoHap2

# Scaffold hap2 contigs
ragtag.py scaffold -u bSpiTri1.HiC.hap2.20240426.fasta legoContigs_hap2.fasta -o ragtag_legoHap2_on_LeGo

# terminal output: 7_scaffoldNuclearContigs.txt
