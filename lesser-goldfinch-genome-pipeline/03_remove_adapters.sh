# Step 2: Remove reads that still have adapters (cutadapt v5.0)

cutadapt \
  -a First_adapter=ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT \
  -a Second_adapter=ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT \
  --error-rate=0.1 \
  --overlap=35 \
  --revcomp \
  --discard-trimmed \
  -o trimmedlego.fastq \
  rawlego.fastq


