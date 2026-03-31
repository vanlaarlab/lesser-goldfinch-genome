# Step 10: Checking for adaptor and foreign DNA contamination (scaffold-level) (FCS)

# Adaptor contamination screen
./run_fcsadaptor.sh --fasta-input ragtag.scaffold.Hap1.fasta --output-dir ./outputdirHap1 --euk
./run_fcsadaptor.sh --fasta-input ragtag.scaffold.Hap2.fasta --output-dir ./outputdirHap2 --euk

# Foreign DNA contamination screen
python3 ./fcs.py screen genome --fasta ragtag.scaffold.lego.Hap1.fasta --out-dir ./gx_outLeGoHap1/ --gx-db "$PWD/db" --tax-id 54559
python3 ./fcs.py screen genome --fasta ragtag.scaffold.Hap2.fasta --out-dir ./gx_outHap2/ --gx-db "$PWD/db" --tax-id 2201931


