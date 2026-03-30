# Step 6a: Assessing nuclear genome completeness (contig-level) (busco v5.8.3)

busco -c 4 -m genome -i legoContigs_hap1.fasta -o buscolegoContigsHap1 -l passeriformes_odb12
busco -c 4 -m genome -i legoContigs_hap2.fasta -o buscolegoContigsHap2 -l passeriformes_odb12

# terminal output: 6_buscolegoContigs.txt
# Step 6b: Visualize quality metrics (blobtoolkit v4.4.6)

docker run --rm -it \
  -v ~/Lego_assembly/contigAssembly/assemblies:/data \
  genomehubs/blobtoolkit:latest \
  bash

cd /data

blobtools create \
  --fasta /data/legoContigs_hap1.fasta \
  /data/blobdir

blobtools add \
  --busco /data/full_table.tsv \
  /data/blobdir

# To host the viewer (expose ports):
docker run --rm -it \
  -p 8001:8001 \
  -p 8002:8002 \
  -v ~/Lego_assembly/contigAssembly/assemblies:/data \
  genomehubs/blobtoolkit:latest \
  bash

blobtools host --api-port 8001 --port 8002 --hostname 0.0.0.0 /data

# View at: http://localhost:8002
