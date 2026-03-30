# Step 9: Assessing nuclear genome completeness (scaffold-level) (busco v5.8.3)

busco -c 4 -m genome -i ragtag.scaffold.Hap1.fasta -o buscolegoScaffoldsHap1 -l passeriformes_odb12
busco -c 4 -m genome -i ragtag.scaffold.Hap2.fasta -o buscolegoScaffoldsHap2 -l passeriformes_odb12

# terminal output: 9_buscolegoScaffolds.txt
