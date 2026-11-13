# Supporting files: Genome Assembly and Annotation of the Lesser Goldfinch (*Spinus psaltria*)

These are the output files behind the numbers reported in the manuscript, numbered in the order the paper uses
them. Each file is the original program output, renamed only. Large files are gzip-compressed (`.gz`).
`CHECKSUMS.md5` lists an MD5 checksum for every file (check with `md5sum -c CHECKSUMS.md5`).

Assembly: GenBank GCA_060154605.1 (primary; BioProject PRJNA1281197) and the alternate assembly
(BioProject PRJNA1281196). Code: https://github.com/vanlaarlab/lesser-goldfinch-assembly.

## Reads

| File | Contents | Supports | Tool |
|---|---|---|---|
| `01_reads_bam_to_fastq_samtools.txt` | conversion of the HiFi BAM to FASTQ | 2,939,006 reads converted | SAMtools 1.21 |
| `02_reads_adapter_trimming_cutadapt.txt` | adapter-filtering report (error rate 0.1, minimum overlap 35 bp, `--revcomp`, `--discard-trimmed`) | only 2 reads contained adapters and were removed | Cutadapt 5.0 |
| `03_reads_raw_read_statistics_seqkit.txt` | read count, total bases, length and quality statistics of the raw reads | 2.94 million reads, 48.8 Gb, mean 16.6 kb, read N50 17.9 kb, 97.8% ≥Q20, 94.7% ≥Q30; coverage 34× and 35× | SeqKit 2.1.0 |

## Assembly

| File | Contents | Supports | Tool |
|---|---|---|---|
| `04_assembly_hifiasm_log.txt` | assembly log, command and run statistics (default parameters) | assembly into two pseudohaplotypes; homozygous coverage peak 38× | hifiasm 0.25.0-r726 |
| `05_assembly_contiguity_quast.txt` | contiguity statistics, primary and alternate assemblies | Table 1: lengths, contigs, largest contig, N50/N90, L50/L90, GC, no Ns | QUAST 5.2.0 |
| `06_assembly_contiguity_quast.tsv` | the same statistics as a tab-separated table | Table 1 | QUAST 5.2.0 |
| `07_assembly_busco_primary.txt` | BUSCO summary, primary assembly (genome mode, miniprot) | 97.8% complete (97.5% single-copy, 0.2% duplicated, 0.1% fragmented, 2.1% missing) | BUSCO 6.1.0, passeriformes_odb12.2 (9,140 genes) |
| `08_assembly_busco_alternate.txt` | BUSCO summary, alternate assembly | 95.5% complete | BUSCO 6.1.0, passeriformes_odb12.2 |

## Sex chromosomes

| File | Contents | Supports | Tool |
|---|---|---|---|
| `09_sexchr_scaffold_assignment.tsv` | every primary scaffold with its alignment to house finch chromosomes and its class (autosomal, Z, W, unplaced, mixed, unaligned) | assignment of 27 scaffolds to the sex chromosomes (22 Z, 5 W) | minimap2 2.24 (`-x asm20`) against the house finch assembly (bHaeMex1) |
| `10_sexchr_assignment_summary.txt` | counts and lengths per class | 22 Z scaffolds (69.6 Mb), 5 W (9.7 Mb), 136 unassigned (9.6 Mb) | custom script |
| `11_sexchr_removed_scaffolds.txt` | the 27 scaffolds excluded as Z or W | repeat analyses without the sex chromosomes | custom script |
| `12_sexchr_read_depth_samtools_coverage.tsv` | per-scaffold read depth from HiFi reads mapped to the primary assembly | Z 19.1×, W 18.0×, autosomes 34.8× (length-weighted means) | minimap2 2.31 (`-ax map-hifi`), SAMtools 1.24 `coverage` |

## Repeats

| File | Contents | Supports | Tool |
|---|---|---|---|
| `13_repeats_curated_library.fa` | curated repeat library used for masking (368 classified families) | repeat annotation; same file as deposited on Dryad | RepeatModeler2 + MCHelper |
| `14_repeats_library_provenance.tsv` | source URL, checksum and family count of the library | library identity | n/a |
| `15_repeats_library_curation_mchelper_log.txt` | MCHelper run on the contigs (`-r A -a F -e 1000`, passeriformes odb12 models) | library curation | MCHelper |
| `16_repeats_repeatmasker_hits.out.gz` | every repeat hit in the primary assembly | basis of all repeat summaries | RepeatMasker 4.1.5 |
| `17_repeats_repeatmasker_summary.tbl` | RepeatMasker's own summary | 31.71% of 1,433,899,291 bp masked (the per-class rows do not read the library's Wicker labels; use file 18) | RepeatMasker 4.1.5 |
| `18_repeats_content_by_chromosome_category.tsv` | repeat content by class for autosomal, Z, W, unassigned, without Z/W and full assembly | Table S2; Z 22.48%, W 80.84%, autosomes 31.68% | custom script |
| `19_repeats_content_with_and_without_ZW.tsv` | each class with and without the sex chromosomes, and its share on Z and W | 31.71% vs 31.84% total; LTR 22.41% vs 22.50% | custom script |
| `20_repeats_divergence_summary_full_assembly.divsum` | Kimura divergence summary, full assembly | landscape source | calcDivergenceFromAlign.pl (RepeatMasker 4.1.5) |
| `21_repeats_divergence_summary_without_ZW.divsum` | Kimura divergence summary, without Z/W | landscape source | calcDivergenceFromAlign.pl |
| `22_repeats_landscape_full_assembly.tsv` | bases per repeat class in 1% Kimura divergence bins, full assembly | landscape comparison | calcDivergenceFromAlign.pl |
| `23_repeats_landscape_without_ZW.tsv` | the same, without Z/W | Figure 1; peak in the 1–2% bin; about 18% of the assembly within 2%; dates in Table 2 | calcDivergenceFromAlign.pl |

In files 22 and 23, each bin is labelled by its lower edge (the script truncates divergence to whole percents), so
bin 1 holds copies with 1–2% divergence.

## Gene annotation

| File | Contents | Supports | Tool |
|---|---|---|---|
| `24_annotation_helixer_gene_models.gff3.gz` | ab initio gene models, vertebrate lineage | 23,684 gene models | Helixer 0.3.7 |
| `25_annotation_gene_model_statistics_agat.txt` | gene model statistics | one transcript per gene; 11.9 exons per transcript; mean gene length 24,720 bp; mean CDS 1,568 bp; 1,054 single-exon genes | AGAT 1.2.0 |
| `26_annotation_predicted_proteins.faa.gz` | predicted protein sequences | 23,684 proteins | AGAT 1.2.0 |
| `27_annotation_busco_proteins.txt` | BUSCO summary, predicted proteome (protein mode) | 94.1% complete (93.8% single-copy, 0.2% duplicated, 1.6% fragmented, 4.3% missing) | BUSCO 6.1.0, passeriformes_odb12.2 |

## Toll-like receptors

Species codes in these files: `hofi` house finch, `pisi` pine siskin, `amgo` American goldfinch, `lego` lesser
goldfinch, `tagu` zebra finch; `chicken` is the reference.

| File | Contents | Supports | Tool |
|---|---|---|---|
| `28_tlr_all_hits.tsv` | every chicken TLR protein alignment in each genome, with identity, similarity, coverage, frameshifts and stops | basis of files 29–33 | miniprot 0.18-r281 |
| `29_tlr_identity_and_similarity.tsv` | best-locus identity and similarity to chicken per TLR and species | identity ranges in the Results | custom script |
| `30_tlr_identity_by_species.tsv` | identity (%) as a TLR × species table | TLR7 79.3–81.2%, TLR3 79.3–79.8%, TLR2B 78.9–79.6%, TLR15 66.2–69.3%; TLR5 71.0–72.2% in finches vs 79.4% in zebra finch | custom script |
| `31_tlr_similarity_by_species.tsv` | similarity (%) as a TLR × species table | Results | custom script |
| `32_tlr_coverage_by_species.tsv` | coverage of the chicken protein (%) | full-length criterion (≥ 50%) | custom script |
| `33_tlr_copy_number_by_species.tsv` | full-length copies / all loci per TLR and species | Table 3 | custom script |
| `34_tlr_loci_and_flanking_genes_all_genomes.tsv` | every TLR locus with its neighbouring genes | gene-neighbourhood comparisons (TLR15, TLR21, TLR5) | custom script (`tlr_flanks.py`) |
| `35_tlr5_loci.tsv` | position of the TLR5 locus in each finch | TLR5 results | miniprot 0.18-r281 |
| `36_tlr5_protein_level_disruptions.tsv` | frameshifts and in-frame stops per species against chicken TLR5 | 8 stops in each finch; 6 frameshifts (American goldfinch) or 7 (others) | miniprot 0.18-r281 |
| `37_tlr5_shared_premature_stop_codons.txt` | stop codons by position, in chicken coordinates | 7 stops shared by all four finches at residues 76, 352, 535, 615, 634, 821 and 829 | MAFFT 7.525, custom script |
| `38_tlr5_codon_alignment_nucleotide.fasta` | codon-aware TLR5 alignment, nucleotides | Figure 2C | MACSE 2.07 |
| `39_tlr5_codon_alignment_protein.fasta` | the same alignment, translated | Figure 2C | MACSE 2.07 |
| `40_tlr5_codon_level_shared_disruptions.txt` | disruptions shared by all four finches, by alignment column | 11 shared disruptions: 7 stops and 4 frameshifts | MACSE 2.07, custom script |
| `41_tlr5_gene_order.tsv` | genes around TLR5 in each species | Figure 2B | miniprot 0.18-r281, NCBI annotations |
| `42_tlr5_exonerate_house_finch.txt` to `45_tlr5_exonerate_lesser_goldfinch.txt` | independent alignment of chicken TLR5 to each finch locus | independent check of the disruptions | Exonerate 2.4.0 (protein2genome) |
| `46_tlr1_cluster_gene_order.tsv` | genes around the TLR1-family cluster in human, chicken and the five passerines | Figure S1; shared flanking genes TBC1D1, KLF3, FAM114A1, TMEM156, KLHL5 | miniprot 0.18-r281, NCBI annotations |
| `47_tlr_pipeline_versions_and_checksums.txt` | tool versions and checksums of inputs and outputs from the TLR pipeline run | software versions (Table S1) | pipeline record |

## Not included

- the RepeatMasker alignment file (15 GB) and soft-masked assembly, which can be regenerated from file 16 and the
  scripts on GitHub;
- duplicate copies (identical files under other names), LaTeX versions of the QUAST report, and macOS system files;
- an American goldfinch QUAST report from a separate analysis;
- the FCS contamination screen output.
