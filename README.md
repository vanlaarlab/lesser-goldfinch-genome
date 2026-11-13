# Lesser Goldfinch (*Spinus psaltria*) Reference Genome

This repository contains the code for *Genome Assembly and Annotation of the Lesser Goldfinch (Spinus psaltria)*
(Cortez, Vargas, Slade and Van Laar). Citation and DOI will be added on publication.

## The study

The lesser goldfinch is a widespread, variably plumaged finch of western North America and northern South
America. The only previous genome assembly for the species was a fragmented short-read assembly. We generated a
high-quality reference genome for a female lesser goldfinch from PacBio HiFi long reads and used it to ask two
questions:

1. **Is the Toll-like receptor gene *TLR5* disrupted, and is the disruption specific to *Spinus*?** TLR5 is the
   main receptor for bacterial flagellin. Lesser and American goldfinches kill *Salmonella* less effectively than
   house finches, so we tested whether *TLR5* loss is restricted to *Spinus* or shared with the house finch.
2. **How much of the genome is repetitive, when were transposable elements (TEs) last active, and how do the sex
   chromosomes affect these estimates?**

Main results: the primary assembly is 1.43 Gb in 938 contigs (N50 12.9 Mb, 97.8% complete BUSCOs). Repeats make up
31.7% of it, mostly LTR retrotransposons; Z-linked sequence is less repetitive than the autosomes and W-linked
sequence is 81% repetitive. The most recent burst of TE activity dates to about 2.0-6.1 million years ago. *TLR5*
carries 11 disruptions shared by the lesser goldfinch, American goldfinch, pine siskin and house finch, so it was
pseudogenized before these species diverged and cannot alone explain the lower *Salmonella* resistance of
*Spinus*.

## How we did it

**Sequencing and assembly.** DNA from a blood sample was sequenced with PacBio HiFi. Reads were converted from
BAM to FASTQ, summarized with SeqKit, and screened for leftover adapters with Cutadapt (2 of 2.94 million reads
were removed). hifiasm assembled two pseudohaplotypes, used as the primary and alternate assemblies, and both were
screened for contamination with NCBI FCS. We assessed contiguity with QUAST and completeness with BUSCO.

**Gene annotation.** Helixer predicted genes ab initio on the primary assembly, without transcript evidence.
AGAT summarized the gene models and extracted their proteins, whose completeness we assessed with BUSCO. The
annotation is exploratory; the TLR analyses below do not depend on it.

**Repeats.** RepeatModeler2 identified repeat families de novo, and MCHelper curated them into a library of 368
classified families. RepeatMasker soft-masked the primary assembly with this library. To identify the sex
chromosomes, we aligned the lesser goldfinch contigs to the house finch assembly, whose Z and W are assembled,
and assigned contigs to Z or W when at least half of their aligned length fell on those chromosomes. Read depth
confirmed the assignment: Z- and W-linked contigs have about half the autosomal depth, as expected in a female.
We summarized repeat content for autosomal, Z and W contigs and for the assembly with and without them. Kimura
divergence of each repeat copy from its family consensus gave the repeat landscape, and the most recent peak was
dated with K = 2μt using pedigree-based mutation rates from three songbirds (collared flycatcher, zebra finch and
great reed warbler), converted with the lesser goldfinch generation length (2.8 years).

**Toll-like receptors.** We aligned chicken TLR proteins to the primary assemblies of the lesser goldfinch,
American goldfinch, pine siskin and house finch, and to the zebra finch, with miniprot. This gave the TLR
repertoire, copy numbers and identity to chicken for each species, and the genes flanking each locus. For *TLR5*,
we identified frameshifts and premature stop codons in protein alignments (miniprot, MAFFT) and codon-aware
alignments (MACSE), checked them with Exonerate, and mapped the disruptions onto the species tree by parsimony to
infer when they arose. We also compared gene order around the TLR1-family cluster in birds with the human
TLR10-TLR1-TLR6 cluster.

## Scripts

The scripts in `scripts/` are numbered in the order they were run:

| Step | Script | Purpose |
|---|---|---|
| 1 | `01_download_raw_reads.sh` | download the HiFi reads from the sequencing core |
| 2 | `02_bam_to_fastq.sh` | convert BAM to FASTQ (SAMtools 1.21) |
| 3 | `03_read_statistics.sh` | read statistics (SeqKit 2.1.0) |
| 4 | `04_remove_adapters.sh` | remove reads containing adapters (Cutadapt 5.0) |
| 5 | `05_assemble_contigs.sh` | assemble the primary and alternate pseudohaplotypes (hifiasm 0.25.0) |
| 6 | `06_contamination_screen.sh` | contamination screen (NCBI FCS 0.5.5) |
| 7 | `07_download_assemblies.sh` | download the GenBank assemblies and reference genomes (NCBI datasets 18.28.0) |
| 8 | `08_quast_assemblies.sh` | contiguity statistics (QUAST 5.2.0) |
| 9 | `09_annotate_genome.sh` | gene prediction (Helixer 0.3.7) |
| 10 | `10_evaluate_annotations.sh` | proteins and gene-model statistics (AGAT 1.2.0) |
| 11 | `11_busco.sh` | completeness of the assemblies and proteins (BUSCO 6.1.0, passeriformes_odb12.2) |
| 12 | `12_repeat_modeling.sh` | de novo repeat families (RepeatModeler 2.0.6) |
| 13 | `13_repeat_curation_mchelper.sh` | curate the repeat library (MCHelper 1.7.1) |
| 14 | `14_repeat_mask_sex_chromosomes.sh` | assign Z/W contigs (minimap2 2.24); mask repeats (RepeatMasker 4.1.5); repeat content with and without Z/W |
| 15 | `15_read_depth.sh` | read depth per contig, checking the Z/W assignment (minimap2 2.31, SAMtools 1.24) |
| 16 | `16_kimura_divergence.sh` | repeat divergence (calcDivergenceFromAlign.pl) |
| 17 | `17_repeat_landscape.R` | repeat landscape figure (R 4.4.3) |
| 18 | `18_te_dating.R` | age of the most recent TE peak (R 4.4.3) |
| 19 | `19_tlr_analysis.sh` | TLR analysis and figures (miniprot 0.18-r281, MAFFT 7.525, Exonerate 2.4.0, MACSE 2.07, AGAT 1.2.0, Python 3.13.11, R) |

Each script's header lists its inputs and outputs, the supporting files it produced (see Data), and, for longer
steps, an example `sbatch` command.

## Running the pipeline

Set the project directory once and run the scripts in order from it:

```bash
export BASE=$PWD
mkdir -p logs
```

Steps 1-6 work in `$BASE/assembly/`. Step 7 downloads the public assemblies the later steps use, so it must run
before steps 8-19. Other dependencies: 13 needs 12; 14 needs 7 and 13; 15 needs 4 and 14; 16 needs 14; 17 and 18
need 14 and 16; 19 needs 7.

The analyses ran on a Linux/SLURM cluster (QUAST, AGAT and RepeatModeler as environment modules; the other tools
in conda environments) and, for the earliest steps and MCHelper, on a workstation with Docker. Partitions, memory
and time in the `sbatch` examples should be adjusted for your system.

Two details for reproducing the published numbers exactly: RepeatModeler samples the genome at random, so steps
12-13 give a similar but not identical library, and the library used in the paper (on Dryad) should be placed at
`te_lego/library/lego_curated_sequences_NR.fa` before step 14. The preparation of the BUSCO profile file for
MCHelper in step 13 was not recorded in the original log and follows the MCHelper documentation.

## Data

| Data | Accession |
|---|---|
| Raw HiFi reads | NCBI BioSample SAMN49551328 |
| Primary assembly | GenBank GCA_060154605.1 (BioProject PRJNA1281197) |
| Alternate assembly | BioProject PRJNA1281196 |
| Curated repeat library | Dryad (link to be added on publication) |

The supporting files (the output files behind each number in the paper, numbered 01-47) are archived with the data
on Dryad, with a README listing what each file shows.

## License

To be added.
