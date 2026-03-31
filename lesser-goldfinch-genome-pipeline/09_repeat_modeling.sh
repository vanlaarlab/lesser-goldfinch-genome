# Step 11: Identify repeat families (repeatmodeler v2.0.6)

# Hap1
BuildDatabase -name legoRepeats ragtag.scaffold.Hap1.fasta
RepeatModeler -database legoRepeats -threads 12 -LTRStruct >& run.out &

# Hap2
BuildDatabase -name legoRepeats legoScaffoldsHap2.trimmed_scafs.fa
RepeatModeler -database legoRepeats -threads 12 -LTRStruct >& run.out &


