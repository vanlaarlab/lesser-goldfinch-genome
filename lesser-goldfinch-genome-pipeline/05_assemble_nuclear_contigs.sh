# Step 4: Assemble nuclear genome (contig-level) (hifiasm v0.25.0, seqtk v1.4)

hifiasm -o legoContigs -t 16 trimmedlego.fastq
awk '/^S/{print ">"$2"\n"$3}' legoContigs.bp.hap1.p_ctg.gfa | fold > legoContigs_hap1.fasta
awk '/^S/{print ">"$2"\n"$3}' legoContigs.bp.hap2.p_ctg.gfa | fold > legoContigs_hap2.fasta

# terminal output: 4_assembleNuclearContigs.txt
