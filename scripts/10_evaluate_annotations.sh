#!/usr/bin/env bash
##################################################################################################################################
# Step 10: Extract predicted proteins and annotation statistics (AGAT v1.2.0)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   25_annotation_gene_model_statistics_agat.txt  gene-model statistics
#   26_annotation_predicted_proteins.faa.gz       predicted proteins
# Input:  Helixer GFF3 and primary assembly from step 9
# Output: eggnog/proteins/lego_pri.faa, helixer/output/lego/lego_pri.stats.txt
# BUSCO on the predicted proteins is run in step 11.
BASE=${BASE:?set BASE to the project directory}
cd "$BASE"
module load agat/1.2.0
mkdir -p eggnog/proteins

for sp in lego; do
  agat_sp_extract_sequences.pl \
    -g helixer/output/$sp/${sp}_pri.helixer.gff3 \
    -f genomes/$sp/${sp}_pri.fna \
    -p \
    -o eggnog/proteins/${sp}_pri.faa
done

ls -lh eggnog/proteins/
grep -c ">" eggnog/proteins/*.faa


agat_sp_statistics.pl \
  --gff helixer/output/lego/lego_pri.helixer.gff3 \
  -o helixer/output/lego/lego_pri.stats.txt