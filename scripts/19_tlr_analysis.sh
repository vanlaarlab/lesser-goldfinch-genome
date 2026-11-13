#!/usr/bin/env bash
##################################################################################################################################
# Step 19: Toll-like receptor analysis
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   28_tlr_all_hits.tsv                                                          every chicken TLR alignment in each genome
#   29_tlr_identity_and_similarity.tsv                                           best-locus identity and similarity
#   30_tlr_identity_by_species.tsv                                               identity by TLR and species
#   31_tlr_similarity_by_species.tsv                                             similarity by TLR and species
#   32_tlr_coverage_by_species.tsv                                               coverage of the chicken protein
#   33_tlr_copy_number_by_species.tsv                                            copy number (Table 3)
#   34_tlr_loci_and_flanking_genes_all_genomes.tsv                               TLR loci and neighboring genes
#   35_tlr5_loci.tsv                                                             TLR5 locus in each finch
#   36_tlr5_protein_level_disruptions.tsv                                        TLR5 frameshifts and stops per species
#   37_tlr5_shared_premature_stop_codons.txt                                     the 7 shared stops in chicken coordinates
#   38_tlr5_codon_alignment_nucleotide.fasta                                     TLR5 codon alignment, nucleotides (Figure 2C)
#   39_tlr5_codon_alignment_protein.fasta                                        TLR5 codon alignment, translated (Figure 2C)
#   40_tlr5_codon_level_shared_disruptions.txt                                   the 11 shared disruptions
#   41_tlr5_gene_order.tsv                                                       genes around TLR5 (Figure 2B)
#   42_tlr5_exonerate_house_finch.txt to 45_tlr5_exonerate_lesser_goldfinch.txt  independent Exonerate alignments
#   46_tlr1_cluster_gene_order.tsv                                               TLR1-family gene order (Figure S1)
#   47_tlr_pipeline_versions_and_checksums.txt                                   tool versions and checksums of the run
# TLR analysis of four finches (miniprot v0.18-r281, MAFFT v7.525, SAMtools v1.24, Exonerate v2.4.0, MACSE v2.07, AGAT v1.2.0, Python v3.13.11)

# Identifies the Toll-like receptor (TLR) genes of the American goldfinch,
# house finch, lesser goldfinch and pine siskin primary assemblies (zebra finch
# as passerine reference) from chicken proteins aligned genome-wide; builds the
# TLR complement table; characterizes TLR5 (protein- and codon-level
# alignments, shared premature stops and frameshifts, exonerate check); and
# compares the TLR1-family cluster with human (T2T-CHM13) and chicken.
#
# Sections
#    1  configuration
#    2  analysis scripts (written to $BASE/scripts, then run by the steps below)
#    3  check       inputs, scripts and conda environments
#    4  refs        chicken CDS into refs/chicken/ (once)
#    5  proteome    chicken proteins aligned genome-wide with miniprot (reused if present)
#    6  flanks      TLR loci and neighboring genes in every genome
#    7  complement  TLR identity, similarity, coverage and copy number
#    8  tlr5        TLR5: locus, chicken-anchored translations, shared stops, exonerate
#    9  codon       TLR5 codon alignment (MACSE), shared frameshifts and stops
#   10  figures     TLR5 figure inputs in figures/tlr5/
#   11  tlr1        TLR1-family cluster synteny in figures/tlr1/
#   12  draw        R figures (here if R has ggplot2, patchwork, gggenes)
#   13  record      tool versions, conda environments, checksums
#   14  summary
#
# Usage, from the project directory:
#   mkdir -p logs
#   BASE=$PWD sbatch --cpus-per-task=16 --mem=48G --time=2-00:00:00 scripts/19_tlr_analysis.sh   # all sections
#   BASE=$PWD bash   scripts/19_tlr_analysis.sh tlr5 codon                                         # selected steps
#   FRESH=1 BASE=$PWD sbatch ... scripts/19_tlr_analysis.sh                                        # remove previous results first
# Inputs below are downloaded by scripts/07_download_assemblies.sh.
#
# Expected under $BASE:
#   genomes/{amgo,hofi,lego,pisi}/<sp>_pri.fna, genomes/tagu/bTaeGut7_mat.fna
#   genomes/tagu/ncbi_dataset/data/GCF_048771995.1/{protein.faa,cds_from_genomic.fna,genomic.gff}
#   gaga_prot.faa, gaga_genomic.gff        chicken GCF_016700215.2 proteins and annotation
#   refs/chicken/cds_from_genomic.fna      chicken CDS (or set CHICKEN_CDS_SRC to copy it)
#   human T2T-CHM13 annotation (GCF_009914755.1): set HUMAN_LOCAL, or it is downloaded
# Requirements: conda environments miniprot (miniprot, samtools, mafft),
#   exonerate_env (exonerate), macse (MACSE v2); AGAT 1.2.0 (module agat/1.2.0);
#   python3; NCBI datasets (only if the human annotation must be downloaded);
#   R with ggplot2, patchwork and gggenes for the figures.
# =============================================================================

# =============================================================================
# Section 1: configuration
# =============================================================================
set -euo pipefail

BASE=${BASE:?set BASE to the project directory}
S=$BASE/scripts
THREADS=${SLURM_CPUS_PER_TASK:-8}
CHICKEN_CDS_SRC=${CHICKEN_CDS_SRC:-}   # chicken GCF_016700215.2 cds_from_genomic.fna; only read if refs/chicken/ is empty
TAGU=genomes/tagu/ncbi_dataset/data/GCF_048771995.1
GENOMES=(genomes/amgo/amgo_pri.fna genomes/hofi/hofi_pri.fna genomes/lego/lego_pri.fna
         genomes/pisi/pisi_pri.fna genomes/tagu/bTaeGut7_mat.fna)
ALL_STEPS="check refs proteome flanks complement tlr5 codon figures tlr1 draw record"
STEPS=${*:-all}
REC=results/tlr_pipeline

die()  { echo "ERROR: $*" >&2; exit 1; }
log()  { echo; echo "[$(date '+%F %T')] ==== $*"; }
want() { [ "$STEPS" = "all" ] || [[ " $STEPS " == *" $1 "* ]]; }

for s in $STEPS; do
  [ "$s" = "all" ] || [[ " $ALL_STEPS " == *" $s "* ]] || die "unknown step '$s' (steps: $ALL_STEPS)"
done

cd "$BASE"
mkdir -p logs "$REC"
module load miniconda3/3.13
eval "$(conda shell.bash hook)"
echo "TLR pipeline on $(hostname), $(date), steps: $STEPS, threads: $THREADS"

if [ "${FRESH:-0}" = "1" ]; then
  log "FRESH=1: removing previous TLR results (miniprot GFFs are kept)"
  rm -rf tlr5 macse tlr_table annotation/tlr_synteny figures/tlr5 figures/tlr1 tlr1_synteny
fi


# =============================================================================
# Section 2: analysis scripts
#   Each script below is written to $BASE/scripts and run by a later section.
# =============================================================================
mkdir -p "$S"

# ---- tlr_flanks.py
cat > "$S/tlr_flanks.py" << '__END_TLR_FLANKS_PY__'
#!/usr/bin/env python3
"""Usage: tlr_flanks.py miniprot.gff chicken.gff sample [flank=2] [min_positive=0.5]"""
import sys, re

mp_gff, chick_gff, sample = sys.argv[1:4]
FLANK  = int(sys.argv[4]) if len(sys.argv) > 4 else 2
MINPOS = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5

def attrs(s):
    return dict(kv.split('=', 1) for kv in s.split(';') if '=' in kv)

p2g = {}
with open(chick_gff) as f:
    for l in f:
        if l.startswith('#'): continue
        c = l.rstrip('\n').split('\t')
        if len(c) < 9 or c[2] != 'CDS': continue
        a = attrs(c[8])
        if 'protein_id' in a and 'gene' in a:
            p2g[a['protein_id']] = a['gene']
p2g = {k: ('TLR2B' if v == 'TLR2' else v) for k, v in p2g.items()}

hits = []
with open(mp_gff) as f:
    for l in f:
        if l.startswith('#'): continue
        c = l.rstrip('\n').split('\t')
        if len(c) < 9 or c[2] != 'mRNA': continue
        a = attrs(c[8])
        if float(a.get('Positive', 0)) < MINPOS: continue
        prot = a['Target'].split()[0]
        hits.append(dict(chr=c[0], s=int(c[3]), e=int(c[4]), strand=c[6],
                         score=float(c[5]) if c[5] != '.' else 0.0,
                         prot=prot, gene=p2g.get(prot, prot),
                         ident=float(a.get('Identity', 0))))

hits.sort(key=lambda h: (h['chr'], h['strand'], h['s']))
loci, cur = [], None
for h in hits:
    if cur and h['chr'] == cur['chr'] and h['strand'] == cur['strand'] and h['s'] <= cur['e']:
        cur['e'] = max(cur['e'], h['e'])
        if h['score'] > cur['best']['score']: cur['best'] = h
    else:
        if cur: loci.append(cur)
        cur = dict(chr=h['chr'], strand=h['strand'], s=h['s'], e=h['e'], best=h)
if cur: loci.append(cur)
loci.sort(key=lambda x: (x['chr'], x['s']))

is_tlr = lambda L: re.match(r'^TLR\d', L['best']['gene']) is not None

print('\t'.join(['sample','tlr_locus','offset','chr','start','end','strand',
                 'gene','chicken_protein','identity','is_tlr']))
n = 0
for i, L in enumerate(loci):
    if not is_tlr(L): continue
    n += 1
    tlr_id = f"{sample}_TLR{n:02d}_{L['best']['gene']}"
    for off in range(-FLANK, FLANK + 1):
        j = i + off
        if j < 0 or j >= len(loci) or loci[j]['chr'] != L['chr']: continue
        M = loci[j]; b = M['best']
        print('\t'.join(map(str, [sample, tlr_id, off, M['chr'], M['s'], M['e'], M['strand'],
                                   b['gene'], b['prot'], b['ident'], int(is_tlr(M))])))
__END_TLR_FLANKS_PY__

# ---- tlr_similarity_table.py
cat > "$S/tlr_similarity_table.py" << '__END_TLR_SIMILARITY_TABLE_PY__'
#!/usr/bin/env python3
"""
TLR similarity table: chicken TLRs (rows) versus every bird (columns).

For each genome, the miniprot alignments of chicken proteins are merged into
loci (overlapping hits on the same strand); each locus is named after the
chicken gene of its best-scoring hit, and loci named after a TLR are kept.
For every chicken TLR and genome the best locus is reported with:

  identity    fraction of aligned residues that are identical  (miniprot Identity)
  similarity  fraction of aligned residues with a positive BLOSUM score
              (miniprot Positive)
  coverage    aligned part of the chicken protein / chicken protein length
  copies      number of loci named after that TLR, as "full-length/all", where
              full-length means coverage >= --min-cov

Identity and similarity are measured over the aligned residues only, so read
them together with coverage: a short fragment can score highly.

Usage (from finches/):
  python3 tlr_similarity_table.py
  python3 tlr_similarity_table.py --pri-only --outdir tlr_table_pri

Inputs:
  gaga_genomic.gff               chicken annotation (protein_id -> gene symbol)
  gaga_prot.faa                  chicken proteins (lengths)
  annotation/miniprot/*/*.gff    miniprot alignments, one file per genome
Only the Python standard library is needed.
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict

PREFERRED = ["tagu", "hofi_pri", "pisi_pri", "amgo_pri", "lego_pri"]


def attrs(s):
    return dict(kv.split("=", 1) for kv in s.split(";") if "=" in kv)


def tlr_key(g):
    m = re.match(r"^TLR(\d+)(.*)$", g)
    return (int(m.group(1)), m.group(2)) if m else (9999, g)


def read_fasta_lengths(path):
    lens, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            lens[name] = 0
        elif name:
            lens[name] += len(line.strip())
    return lens


def read_chicken(gff):
    """protein_id -> gene symbol, and gene symbol -> chr:start-end."""
    p2g, loc = {}, {}
    for line in open(gff):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9:
            continue
        if c[2] == "CDS":
            a = attrs(c[8])
            if "protein_id" in a and "gene" in a:
                p2g[a["protein_id"]] = a["gene"]
        elif c[2] == "gene":
            a = attrs(c[8])
            g = a.get("gene")
            if g and g not in loc:
                loc[g] = f"{c[0]}:{c[3]}-{c[4]}"
    # NCBI's symbol for the second chicken TLR2 copy is the bare "TLR2"
    for k, v in list(p2g.items()):
        if v == "TLR2":
            p2g[k] = "TLR2B"
    if "TLR2" in loc:
        loc["TLR2B"] = loc.pop("TLR2")
    return p2g, loc


def genome_name(path):
    base = os.path.basename(path)[:-4]
    return "tagu" if base.startswith("bTaeGut7") else base


def load_hits(path, p2g, plen):
    hits = []
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "mRNA":
            continue
        a = attrs(c[8])
        tgt = a["Target"].split()
        prot, qs, qe = tgt[0], int(tgt[1]), int(tgt[2])
        hits.append({
            "chr": c[0], "s": int(c[3]), "e": int(c[4]), "strand": c[6],
            "score": float(c[5]) if c[5] != "." else 0.0,
            "prot": prot, "gene": p2g.get(prot, prot),
            "ident": float(a.get("Identity", 0)),
            "sim": float(a.get("Positive", 0)),
            "cov": (qe - qs + 1) / plen[prot] if prot in plen else float("nan"),
            "fs": int(a.get("Frameshift", 0)),
            "stops": int(a.get("StopCodon", 0)),
        })
    return hits


def tlr_loci(hits):
    """Merge overlapping same-strand hits into loci, keep those named TLR*."""
    hits = sorted(hits, key=lambda h: (h["chr"], h["strand"], h["s"]))
    loci, cur = [], None
    for h in hits:
        if (cur and h["chr"] == cur["chr"] and h["strand"] == cur["strand"]
                and h["s"] <= cur["e"]):
            cur["e"] = max(cur["e"], h["e"])
            if h["score"] > cur["best"]["score"]:
                cur["best"] = h
        else:
            if cur:
                loci.append(cur)
            cur = {"chr": h["chr"], "strand": h["strand"], "s": h["s"],
                   "e": h["e"], "best": h}
    if cur:
        loci.append(cur)
    return [L for L in loci if re.match(r"^TLR\d", L["best"]["gene"])]


def pct(x):
    return "NA" if x != x else f"{100 * x:.1f}"


def write_wide(path, genes, genomes, fn, first_cols=None):
    with open(path, "w") as out:
        head = ["chicken_gene"] + (list(first_cols[1]) if first_cols else []) + genomes
        out.write("\t".join(head) + "\n")
        for g in genes:
            row = [g]
            if first_cols:
                row += [str(first_cols[0][g].get(k, "")) for k in first_cols[1]]
            row += [fn(g, gen) for gen in genomes]
            out.write("\t".join(row) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chicken-gff", default="gaga_genomic.gff")
    ap.add_argument("--chicken-faa", default="gaga_prot.faa")
    ap.add_argument("--miniprot-dir", default="annotation/miniprot")
    ap.add_argument("--outdir", default="tlr_table")
    ap.add_argument("--min-cov", type=float, default=0.5,
                    help="coverage needed to count a locus as a full-length copy")
    ap.add_argument("--pri-only", action="store_true",
                    help="primary assemblies and zebra finch only (no alt haplotypes)")
    args = ap.parse_args()

    for f in (args.chicken_gff, args.chicken_faa):
        if not os.path.isfile(f):
            sys.exit(f"ERROR: not found: {f}")
    files = sorted(glob.glob(os.path.join(args.miniprot_dir, "*", "*.gff")))
    if not files:
        sys.exit(f"ERROR: no GFFs under {args.miniprot_dir}/*/")

    p2g, gene_loc = read_chicken(args.chicken_gff)
    plen = read_fasta_lengths(args.chicken_faa)

    # representative chicken protein per TLR gene = the longest isoform
    rep = {}
    for prot, g in p2g.items():
        if re.match(r"^TLR\d", g) and prot in plen:
            if g not in rep or plen[prot] > plen[rep[g]]:
                rep[g] = prot
    genes = sorted(rep, key=tlr_key)
    print("chicken TLR genes:", ", ".join(genes))

    results = {}                      # genome -> gene -> list of loci
    for f in files:
        gen = genome_name(f)
        if args.pri_only and not (gen == "tagu" or gen.endswith("_pri")):
            continue
        by_gene = defaultdict(list)
        for L in tlr_loci(load_hits(f, p2g, plen)):
            by_gene[L["best"]["gene"]].append(L)
        results[gen] = by_gene
        print(f"{gen:14s} {sum(len(v) for v in by_gene.values())} TLR loci")

    genomes = ([g for g in PREFERRED if g in results]
               + sorted(g for g in results if g not in PREFERRED))

    def best(g, gen):
        loci = results[gen].get(g, [])
        return max(loci, key=lambda L: L["best"]["score"]) if loci else None

    os.makedirs(args.outdir, exist_ok=True)

    # long table: one row per chicken TLR x genome (best locus)
    with open(os.path.join(args.outdir, "tlr_hits_long.tsv"), "w") as out:
        out.write("\t".join([
            "chicken_gene", "chicken_protein", "chicken_len_aa", "genome",
            "contig", "start", "end", "strand", "identity_pct",
            "similarity_pct", "coverage_pct", "frameshifts", "stop_codons",
            "loci_full_length", "loci_total"]) + "\n")
        for g in genes:
            for gen in genomes:
                L = best(g, gen)
                loci = results[gen].get(g, [])
                nfull = sum(1 for x in loci if x["best"]["cov"] >= args.min_cov)
                if L is None:
                    row = [g, rep[g], str(plen[rep[g]]), gen] + ["NA"] * 9 + ["0", "0"]
                else:
                    b = L["best"]
                    row = [g, rep[g], str(plen[rep[g]]), gen, L["chr"],
                           str(L["s"]), str(L["e"]), L["strand"],
                           pct(b["ident"]), pct(b["sim"]), pct(b["cov"]),
                           str(b["fs"]), str(b["stops"]),
                           str(nfull), str(len(loci))]
                out.write("\t".join(row) + "\n")

    info = ({g: {"chicken_protein": rep[g], "chicken_len_aa": plen[rep[g]],
                 "chicken_location": gene_loc.get(g, "")} for g in genes},
            ["chicken_protein", "chicken_len_aa", "chicken_location"])

    def cell(key):
        def f(g, gen):
            L = best(g, gen)
            return "NA" if L is None else pct(L["best"][key])
        return f

    def copies(g, gen):
        loci = results[gen].get(g, [])
        nfull = sum(1 for x in loci if x["best"]["cov"] >= args.min_cov)
        return f"{nfull}/{len(loci)}"

    write_wide(os.path.join(args.outdir, "tlr_identity_wide.tsv"),
               genes, genomes, cell("ident"), info)
    write_wide(os.path.join(args.outdir, "tlr_similarity_wide.tsv"),
               genes, genomes, cell("sim"), info)
    write_wide(os.path.join(args.outdir, "tlr_coverage_wide.tsv"),
               genes, genomes, cell("cov"))
    write_wide(os.path.join(args.outdir, "tlr_copies_wide.tsv"),
               genes, genomes, copies)

    def combined(g, gen):
        L = best(g, gen)
        if L is None:
            return "-"
        b = L["best"]
        return f"{pct(b['ident'])} / {pct(b['sim'])}"
    write_wide(os.path.join(args.outdir, "tlr_identity_similarity.tsv"),
               genes, genomes, combined)

    # readable summary on screen
    print("\nidentity / similarity (%) of the best locus; "
          "copies = full-length/all loci")
    w = max(len(x) for x in genomes + ["chicken_gene"]) + 2
    print("chicken_gene".ljust(14) + "".join(x.ljust(max(w, 16)) for x in genomes))
    for g in genes:
        print(g.ljust(14) + "".join(combined(g, x).ljust(max(w, 16)) for x in genomes))
    print("\ncopies (full-length/all)")
    for g in genes:
        print(g.ljust(14) + "".join(copies(g, x).ljust(max(w, 16)) for x in genomes))
    print(f"\nwrote tables to {args.outdir}/")


if __name__ == "__main__":
    main()
__END_TLR_SIMILARITY_TABLE_PY__

# ---- tlr5_analysis.sh
cat > "$S/tlr5_analysis.sh" << '__END_TLR5_ANALYSIS_SH__'
#!/bin/bash
# =============================================================================
# TLR5 pseudogenization in four fringillid finches
#
# Tests whether TLR5 is intact or pseudogenized in American goldfinch (amgo),
# house finch (hofi), lesser goldfinch (lego) and pine siskin (pisi), using
# chicken and zebra finch TLR5 as intact references.
#
# Logic:
#   1. Locate the TLR5 locus in each genome with a passerine query (zebra finch)
#   2. Align the chicken TLR5 protein to all five passerine genomes, so every
#      translation sits in the same reference frame and is directly comparable
#   3. Align those translations and find premature stop codons shared across
#      species -> shared positions imply one ancestral loss, not independent ones
#   4. Confirm with a second aligner (exonerate), which also reports frameshifts
#
# Inputs (relative to BASE):
#   genomes/<sp>/<sp>_pri.fna                                 finch assemblies
#   genomes/tagu/bTaeGut7_mat.fna                             zebra finch assembly
#   genomes/tagu/ncbi_dataset/data/GCF_048771995.1/protein.faa
#   gaga_prot.faa    chicken GRCg7w proteins (GCF_016700215.2)
#
# Environments: miniprot (miniprot, samtools, mafft), exonerate_env (exonerate)
#
# Runs in a few minutes on the login node. To queue it instead:
#   sbatch -p windfall -c 8 --mem=16G -t 02:00:00 \
#          -o logs/tlr5_%j.out -e logs/tlr5_%j.err \
#          --wrap "bash scripts/tlr5_analysis.sh"
# =============================================================================

set -euo pipefail

BASE=${BASE:?set BASE to the project directory}
TAGU_REF=$BASE/genomes/tagu/ncbi_dataset/data/GCF_048771995.1
SPECIES=(amgo hofi lego pisi)
THREADS=8

die() { echo "ERROR: $*" >&2; exit 1; }

cd "$BASE" || die "cannot cd to $BASE"
mkdir -p tlr5/queries tlr5/loci tlr5/aln tlr5/results

module load miniconda3/3.13
eval "$(conda shell.bash hook)"
conda activate miniprot


# -----------------------------------------------------------------------------
# 1. Reference proteins
#    XP_072783695.1  zebra finch TLR5, 859 aa  (intact passerine, used to locate)
#    NP_001384988.1  chicken TLR5,    881 aa  (common frame for all species)
# -----------------------------------------------------------------------------
echo "[1] extracting reference proteins"

awk '/^>/{p=(index($0,">XP_072783695.1")==1)} p' "$TAGU_REF/protein.faa" \
  > tlr5/queries/tagu_TLR5.faa
awk '/^>/{p=(index($0,">NP_001384988.1")==1)} p' gaga_prot.faa \
  > tlr5/queries/gaga_TLR5.faa

for f in tlr5/queries/tagu_TLR5.faa tlr5/queries/gaga_TLR5.faa; do
  [ "$(grep -c '>' "$f")" -eq 1 ] || die "$f is not a single sequence"
done


# -----------------------------------------------------------------------------
# 2. Locate the TLR5 locus in each finch (zebra finch query, closest intact ref)
#
#    Each GFF is written to a .tmp and renamed, so no later step can read a
#    partially written file.
# -----------------------------------------------------------------------------
echo "[2] locating TLR5 in each finch assembly"

for sp in "${SPECIES[@]}"; do
  g="genomes/$sp/${sp}_pri.fna"
  [ -s "$g" ] || die "genome not found: $g"
  miniprot -t $THREADS --gff --outs=0.5 "$g" tlr5/queries/tagu_TLR5.faa \
    > "tlr5/loci/${sp}_TLR5.gff.tmp" 2> "tlr5/loci/${sp}_TLR5.log" \
    || die "miniprot failed for $sp (see tlr5/loci/${sp}_TLR5.log)"
  mv "tlr5/loci/${sp}_TLR5.gff.tmp" "tlr5/loci/${sp}_TLR5.gff"
done

# locus table; a species with no mRNA line is a hard error, not a silent gap
{
  printf "species\tscaffold\tstart\tend\tstrand\tidentity\tframeshifts\tstops\n"
  for sp in "${SPECIES[@]}"; do
    row=$(awk -F'\t' -v s="$sp" '$3=="mRNA"{
      id=""; fs=""; st=""
      n=split($9,a,";")
      for(i=1;i<=n;i++){
        split(a[i],kv,"=")
        if(kv[1]=="Identity")   id=kv[2]
        if(kv[1]=="Frameshift") fs=kv[2]
        if(kv[1]=="StopCodon")  st=kv[2]
      }
      print s"\t"$1"\t"$4"\t"$5"\t"$7"\t"id"\t"fs"\t"st; exit
    }' "tlr5/loci/${sp}_TLR5.gff")
    [ -n "$row" ] || die "no TLR5 alignment found for $sp"
    echo "$row"
  done
} > tlr5/results/tlr5_loci.tsv

column -t tlr5/results/tlr5_loci.tsv


# -----------------------------------------------------------------------------
# 3. Chicken TLR5 against all five passerine genomes
#    --trans --aln emits a "##STA" line holding the translated CDS of the
#    alignment, with in-frame stops as '*'. Aligning every species to the same
#    chicken query puts all translations in one comparable frame.
# -----------------------------------------------------------------------------
echo "[3] aligning chicken TLR5 to all five passerine genomes"

declare -A GENOME=(
  [amgo]=genomes/amgo/amgo_pri.fna
  [hofi]=genomes/hofi/hofi_pri.fna
  [lego]=genomes/lego/lego_pri.fna
  [pisi]=genomes/pisi/pisi_pri.fna
  [tagu]=genomes/tagu/bTaeGut7_mat.fna          # intact control
)
ALL=(amgo hofi lego pisi tagu)

for sp in "${ALL[@]}"; do
  g=${GENOME[$sp]}
  [ -s "$g" ] || die "genome not found: $g"
  miniprot --trans --aln -t $THREADS "$g" tlr5/queries/gaga_TLR5.faa \
    > "tlr5/aln/${sp}_vs_chicken.txt.tmp" 2> "tlr5/aln/${sp}_vs_chicken.log" \
    || die "miniprot --aln failed for $sp"
  mv "tlr5/aln/${sp}_vs_chicken.txt.tmp" "tlr5/aln/${sp}_vs_chicken.txt"
done

# frameshift / in-frame stop counts from the PAF line (fs:i: and st:i: tags)
{
  printf "species\tframeshifts\tin_frame_stops\n"
  for sp in "${ALL[@]}"; do
    row=$(grep -m1 "^NP_001384988.1" "tlr5/aln/${sp}_vs_chicken.txt" \
      | awk -v s="$sp" '{
          fs=""; st=""
          for(i=1;i<=NF;i++){
            if($i ~ /^fs:i:/) fs=substr($i,6)
            if($i ~ /^st:i:/) st=substr($i,6)
          }
          print s"\t"fs"\t"st
        }')
    [ -n "$row" ] || die "no PAF line for $sp"
    echo "$row"
  done
} > tlr5/results/tlr5_disruptions.tsv

column -t tlr5/results/tlr5_disruptions.tsv


# -----------------------------------------------------------------------------
# 4. Align the translations, plus the chicken query itself as a second control
#    MAFFT drops non-standard characters by default, so --anysymbol is needed
#    to keep the '*' stop codons.
# -----------------------------------------------------------------------------
echo "[4] aligning translations"

for sp in "${ALL[@]}"; do
  seq=$(awk -F'\t' '$1=="##STA"{print $2; exit}' "tlr5/aln/${sp}_vs_chicken.txt")
  [ -n "$seq" ] || die "no ##STA translation for $sp"
  printf ">%s\n%s\n" "$sp" "$seq"
done > tlr5/aln/tlr5_translations.faa

cat tlr5/queries/gaga_TLR5.faa >> tlr5/aln/tlr5_translations.faa

mafft --anysymbol --auto tlr5/aln/tlr5_translations.faa \
  > tlr5/aln/tlr5_translations.mafft.faa 2> tlr5/aln/mafft.log \
  || die "mafft failed (see tlr5/aln/mafft.log)"

grep -q '\*' tlr5/aln/tlr5_translations.mafft.faa \
  || die "no '*' survived alignment; check that --anysymbol is supported"


# -----------------------------------------------------------------------------
# 5. Premature stops: per species, shared across finches, in chicken coordinates
# -----------------------------------------------------------------------------
echo "[5] scoring premature stop codons"

python3 - << 'PYEOF' | tee tlr5/results/tlr5_stops.txt
seqs, name = {}, None
for line in open('tlr5/aln/tlr5_translations.mafft.faa'):
    if line.startswith('>'):
        name = line[1:].split()[0]; seqs[name] = []
    else:
        seqs[name].append(line.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}

finches = ['amgo', 'hofi', 'lego', 'pisi']
chicken = [k for k in seqs if k.startswith('NP_')][0]
L = len(seqs[chicken])

# alignment column -> chicken residue number
col2chick, n = {}, 0
for i, c in enumerate(seqs[chicken]):
    if c != '-':
        n += 1
    col2chick[i + 1] = n

print("stop codons per sequence (alignment columns; the last is the terminal stop)")
for k in finches + ['tagu', chicken]:
    pos = [i + 1 for i, c in enumerate(seqs[k]) if c == '*']
    print(f"  {k:16s} {len(pos):2d}  {pos}")

shared = [i + 1 for i in range(L)
          if all(seqs[f][i] == '*' for f in finches)
          and seqs[chicken][i] not in ('*', '-')]

print(f"\npremature stops shared by all four finches: {len(shared)}")
print(f"{'aln_col':>8}  {'chicken_res':>11}  {'chicken_aa':>10}  {'tagu_aa':>8}")
for col in shared:
    print(f"{col:8d}  {col2chick[col]:11d}  {seqs[chicken][col-1]:>10}  {seqs['tagu'][col-1]:>8}")

print("\nstops not shared by all four (terminal stop included):")
for f in finches:
    p = sorted(set(i + 1 for i, c in enumerate(seqs[f]) if c == '*') - set(shared))
    print(f"  {f}: {p}")
PYEOF


# -----------------------------------------------------------------------------
# 6. Independent confirmation with exonerate
#    '#' marks frameshifts and '***' marks in-frame stops in its alignment.
#    Locus windows come from the step-2 coordinates with 5 kb of flank.
# -----------------------------------------------------------------------------
echo "[6] extracting loci and running exonerate"

while read -r sp scaf start end strand rest; do
  [ "$sp" = "species" ] && continue
  s=$(( start > 5000 ? start - 5000 : 1 ))
  e=$(( end + 5000 ))
  samtools faidx "genomes/$sp/${sp}_pri.fna" "${scaf}:${s}-${e}" \
    > "tlr5/loci/${sp}_locus.fa" || die "samtools faidx failed for $sp"
  [ -s "tlr5/loci/${sp}_locus.fa" ] || die "empty locus file for $sp"
done < tlr5/results/tlr5_loci.tsv

conda activate exonerate_env
for sp in "${SPECIES[@]}"; do
  exonerate --model protein2genome --showalignment yes --showvulgar no \
    tlr5/queries/gaga_TLR5.faa "tlr5/loci/${sp}_locus.fa" \
    > "tlr5/aln/${sp}_exonerate.txt" || die "exonerate failed for $sp"
done

echo
echo "done."
echo "  tlr5/results/tlr5_loci.tsv             locus coordinates per species"
echo "  tlr5/results/tlr5_disruptions.tsv      frameshift and stop counts"
echo "  tlr5/results/tlr5_stops.txt            shared stops in chicken coordinates"
echo "  tlr5/aln/tlr5_translations.mafft.faa   alignment for AliView"
echo "  tlr5/aln/*_exonerate.txt               second-aligner confirmation"
__END_TLR5_ANALYSIS_SH__

# ---- tlr5_macse.sh
cat > "$S/tlr5_macse.sh" << '__END_TLR5_MACSE_SH__'
#!/bin/bash
# =============================================================================
# Codon-aware TLR5 alignment and shared-disruption detection with MACSE v2
#
# Follows the MACSE recommendation for datasets mixing functional genes and
# pseudogenes: functional references go in -seq with default (high) penalties,
# putative pseudogenes go in -seq_lr with lower frameshift/stop penalties, so
# MACSE will place those events rather than mangling the codon structure.
#
#   reliable   (-seq)     chicken + zebra finch TLR5 CDS  (annotated, functional)
#   less rel.  (-seq_lr)  the four finch TLR5 loci        (putative pseudogenes)
#
# Note: alignSequences ignores sequence orientation, so every input must be in
# the forward coding direction. AGAT handles this when extracting from the
# miniprot GFFs (lego is on the minus strand).
#
# Prerequisite: scripts/tlr5_analysis.sh has been run, so tlr5/loci/*_TLR5.gff
# exists.
# =============================================================================

set -euo pipefail

BASE=${BASE:?set BASE to the project directory}
SPECIES=(amgo hofi lego pisi)
GAGA_CDS=$BASE/refs/chicken/cds_from_genomic.fna   # copied from scratch by run_tlr_pipeline.sh
TAGU_CDS=$BASE/genomes/tagu/ncbi_dataset/data/GCF_048771995.1/cds_from_genomic.fna

die() { echo "ERROR: $*" >&2; exit 1; }

cd "$BASE" || die "cannot cd to $BASE"
mkdir -p macse/input macse/out

module load miniconda3/3.13
eval "$(conda shell.bash hook)"


# -----------------------------------------------------------------------------
# 1. Reference CDS (functional, "reliable")
#    NCBI cds_from_genomic.fna headers carry [protein_id=...]
# -----------------------------------------------------------------------------
echo "[1] extracting reference CDS"

[ -s "$GAGA_CDS" ] || die "chicken CDS not found: $GAGA_CDS"
[ -s "$TAGU_CDS" ] || die "zebra finch CDS not found: $TAGU_CDS"

awk '/^>/{p=(index($0,"protein_id=NP_001384988.1")>0); if(p) print ">chicken"; next} p' \
  "$GAGA_CDS" > macse/input/refs.fna
awk '/^>/{p=(index($0,"protein_id=XP_072783695.1")>0); if(p) print ">zebrafinch"; next} p' \
  "$TAGU_CDS" >> macse/input/refs.fna

[ "$(grep -c '>' macse/input/refs.fna)" -eq 2 ] \
  || die "expected 2 reference CDS, got $(grep -c '>' macse/input/refs.fna)"


# -----------------------------------------------------------------------------
# 2. Finch CDS from the miniprot alignments (putative pseudogenes)
#    AGAT concatenates exons in transcript order and reverse-complements
#    minus-strand loci, so the output is forward coding sequence.
# -----------------------------------------------------------------------------
echo "[2] extracting finch CDS"

module load agat/1.2.0
: > macse/input/finches.fna
for sp in "${SPECIES[@]}"; do
  g="tlr5/loci/${sp}_TLR5.gff"
  [ -s "$g" ] || die "missing $g (run scripts/tlr5_analysis.sh first)"
  agat_sp_extract_sequences.pl \
    -g "$g" -f "genomes/$sp/${sp}_pri.fna" -t cds \
    -o "macse/input/${sp}_cds.fna" > "macse/input/${sp}_agat.log" 2>&1 \
    || die "AGAT failed for $sp"
  # single record per locus, renamed to the species
  awk -v s="$sp" 'BEGIN{printf ">%s\n", s} !/^>/{printf "%s", $0} END{print ""}' \
    "macse/input/${sp}_cds.fna" >> macse/input/finches.fna
done

echo "input CDS lengths (expect ~2600 bp for an ~860 aa protein):"
awk '/^>/{if(n)printf "  %-12s %d\n", h, n; h=substr($1,2); n=0; next}{n+=length($0)}
     END{printf "  %-12s %d\n", h, n}' macse/input/refs.fna macse/input/finches.fna


# -----------------------------------------------------------------------------
# 3. MACSE alignment
#    Defaults penalize frameshifts and stops heavily; -fs_lr / -stop_lr lower
#    that for the finch sequences only.
# -----------------------------------------------------------------------------
echo "[3] running MACSE"

conda activate macse
MACSE="macse"                      # conda wrapper; or: java -Xmx4g -jar /path/macse_v2.XX.jar

$MACSE -prog alignSequences \
  -seq    macse/input/refs.fna \
  -seq_lr macse/input/finches.fna \
  -fs_lr 10 -stop_lr 10 \
  -out_NT macse/out/tlr5_NT.fasta \
  -out_AA macse/out/tlr5_AA.fasta \
  > macse/out/macse.log 2>&1 || die "MACSE failed (see macse/out/macse.log)"

[ -s macse/out/tlr5_NT.fasta ] || die "MACSE produced no NT alignment"
grep -c '>' macse/out/tlr5_NT.fasta


# -----------------------------------------------------------------------------
# 4. Score shared frameshifts and premature stops
# -----------------------------------------------------------------------------
echo "[4] scoring disruptions"

python3 scripts/shared_disruptions.py macse/out/tlr5_NT.fasta \
  --refs chicken,zebrafinch \
  | tee macse/out/tlr5_disruptions.txt

echo
echo "done."
echo "  macse/out/tlr5_NT.fasta        codon alignment (view in SeaView or AliView)"
echo "  macse/out/tlr5_AA.fasta        translated alignment"
echo "  macse/out/tlr5_disruptions.txt shared frameshifts and stops"
__END_TLR5_MACSE_SH__

# ---- shared_disruptions.py
cat > "$S/shared_disruptions.py" << '__END_SHARED_DISRUPTIONS_PY__'
#!/usr/bin/env python3
"""
Find frameshifts and premature stop codons in a MACSE nucleotide alignment
and report which sequences share each one.

MACSE marks frameshifts with '!' in its _NT output. Because MACSE preserves
codon structure, homologous codons occupy the same alignment columns, so a
shared stop is an exact column match. Frameshift placement is less certain
(an indel in a homopolymer can be placed at several positions), so only
frameshifts are merged across a small window.

Usage:
  shared_disruptions.py aln_NT.fasta [--refs chicken,zebrafinch] [--fs-window 3]

Any event also present in a reference sequence is flagged, since a "shared"
disruption that appears in a known functional gene is an alignment artifact
rather than a biological one.
"""

import argparse
from collections import defaultdict

STOPS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def find_events(seq):
    """Return [(type, 1-based alignment column), ...] for one aligned sequence."""
    events = []

    # frameshifts: a run of '!' is one event, recorded at its first column
    i = 0
    while i < len(seq):
        if seq[i] == "!":
            start = i
            while i < len(seq) and seq[i] == "!":
                i += 1
            events.append(("frameshift", start + 1))
        else:
            i += 1

    # stop codons, on the codon grid MACSE guarantees
    codons = [(j, seq[j:j + 3]) for j in range(0, len(seq) - 2, 3)]
    real = [(j, c) for j, c in codons if set(c) - set("-")]
    for idx, (j, c) in enumerate(real):
        if c in STOPS:
            # the final codon is the normal terminal stop, not a disruption
            if idx == len(real) - 1:
                continue
            events.append(("stop", j + 1))
    return events


def cluster(flat, fs_window):
    """Group identical-type events across sequences.

    Stops must match exactly. Frameshifts may be merged within fs_window
    columns of the cluster's FIRST position, so clusters cannot chain
    indefinitely across a long run of nearby events.
    """
    window = {"stop": 0, "frameshift": fs_window}
    clusters = []  # [type, first, last, {names}]
    for t, p, n in flat:
        for c in clusters:
            if c[0] == t and p - c[1] <= window[t]:
                c[2] = max(c[2], p)
                c[3].add(n)
                break
        else:
            clusters.append([t, p, p, {n}])
    return clusters


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("alignment")
    ap.add_argument("--refs", default="",
                    help="comma-separated names of known-functional sequences")
    ap.add_argument("--fs-window", type=int, default=3,
                    help="merge frameshifts within this many columns (default 3)")
    args = ap.parse_args()

    seqs = read_fasta(args.alignment)
    refs = {r for r in args.refs.split(",") if r}
    missing = refs - set(seqs)
    if missing:
        raise SystemExit(f"ERROR: reference(s) not in alignment: {sorted(missing)}")
    queries = [n for n in seqs if n not in refs]

    by_seq = {name: find_events(s) for name, s in seqs.items()}
    flat = sorted((t, p, n) for n, ev in by_seq.items() for t, p in ev)
    clusters = cluster(flat, args.fs_window)

    print("All events (a cluster containing a reference is likely an artifact)")
    print("type\tcolumn\tn_seqs\tin_ref\tsequences")
    for t, p1, p2, names in sorted(clusters, key=lambda c: (c[1], c[0])):
        col = str(p1) if p1 == p2 else f"{p1}-{p2}"
        in_ref = "YES" if names & refs else "-"
        print(f"{t}\t{col}\t{len(names)}\t{in_ref}\t{','.join(sorted(names))}")

    if queries:
        shared = [c for c in clusters
                  if set(queries) <= c[3] and not (c[3] & refs)]
        print(f"\nShared by all {len(queries)} query sequences "
              f"({', '.join(sorted(queries))}) and absent from references: "
              f"{len(shared)}")
        print("type\tcolumn")
        for t, p1, p2, _ in sorted(shared, key=lambda c: c[1]):
            print(f"{t}\t{p1 if p1 == p2 else f'{p1}-{p2}'}")

    print("\nPer-sequence summary")
    print("sequence\tframeshifts\tinternal_stops\tnote")
    for name in sorted(seqs):
        ev = by_seq[name]
        fs = sum(1 for t, _ in ev if t == "frameshift")
        st = sum(1 for t, _ in ev if t == "stop")
        note = "reference" if name in refs else (
            "no disruptions detected" if not ev else "disrupted")
        print(f"{name}\t{fs}\t{st}\t{note}")


if __name__ == "__main__":
    main()
__END_SHARED_DISRUPTIONS_PY__

# ---- tlr5_figure_inputs.sh
cat > "$S/tlr5_figure_inputs.sh" << '__END_TLR5_FIGURE_INPUTS_SH__'
#!/bin/bash
# =============================================================================
# Gather inputs for the TLR5 figure. Run on the cluster from finches/.
#
# Produces, in figures/tlr5/:
#   tlr5_synteny.tsv   TLR5 plus two protein-coding neighbors on each side for
#                      chicken, zebra finch and the four finch primary assemblies
#   tlr5_NT.fasta      the MACSE codon alignment (copied from macse/out/)
#
# Neighbor sources:
#   chicken, zebra finch  NCBI annotations (protein-coding genes, gene= symbols)
#   finches               miniprot hits of chicken proteins, from
#                         annotation/tlr_synteny/all_tlr_flanks.tsv (window 2),
#                         labeled with chicken gene symbols
# =============================================================================

set -euo pipefail

BASE=${BASE:?set BASE to the project directory}
cd "$BASE"
OUTDIR=figures/tlr5
mkdir -p "$OUTDIR"

die() { echo "ERROR: $*" >&2; exit 1; }

GAGA_GFF=gaga_genomic.gff
TAGU_GFF=genomes/tagu/ncbi_dataset/data/GCF_048771995.1/genomic.gff
FLANKS=annotation/tlr_synteny/all_tlr_flanks.tsv
MACSE_NT=macse/out/tlr5_NT.fasta

for f in "$GAGA_GFF" "$TAGU_GFF" "$FLANKS" "$MACSE_NT"; do
  [ -s "$f" ] || die "missing input: $f"
done

# TLR5 and two protein-coding genes either side, from an NCBI GFF
flanks_from_gff () {
  local gff=$1 label=$2
  awk -F'\t' '$3=="gene" && $9 ~ /gene_biotype=protein_coding/ {
      n=split($9,a,";"); g=""
      for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]=="gene") g=kv[2]}
      if(g!="") print $1"\t"$4"\t"$5"\t"$7"\t"g
    }' "$gff" \
  | sort -t$'\t' -k1,1 -k2,2n \
  | awk -F'\t' -v lab="$label" '
      {c[NR]=$1; s[NR]=$2; e[NR]=$3; st[NR]=$4; g[NR]=$5}
      END{
        for(i=1;i<=NR;i++) if(g[i]=="TLR5"){
          for(j=i-2;j<=i+2;j++) if(j>=1 && j<=NR && c[j]==c[i])
            print lab"\t"g[j]"\t"c[j]"\t"s[j]"\t"e[j]"\t"st[j]"\t"(j==i?1:0)
        }
      }'
}

{
  printf "species\tgene\tcontig\tstart\tend\tstrand\tfocal\n"
  flanks_from_gff "$GAGA_GFF" chicken
  flanks_from_gff "$TAGU_GFF" zebrafinch
  # flank table columns: sample tlr_locus offset chr start end strand gene ...
  awk -F'\t' 'NR>1 && $1 ~ /_pri$/ && $2 ~ /_TLR5$/ {
      split($1,a,"_")
      print a[1]"\t"$8"\t"$4"\t"$5"\t"$6"\t"$7"\t"($3==0?1:0)
    }' "$FLANKS"
} > "$OUTDIR/tlr5_synteny.tsv"

cp "$MACSE_NT" "$OUTDIR/tlr5_NT.fasta"

echo "genes per species (expect 5, with one TLR5 each):"
awk -F'\t' 'NR>1{n[$1]++; f[$1]+=$7}
  END{for(s in n) printf "  %-12s %d genes, %d TLR5\n", s, n[s], f[s]}' \
  "$OUTDIR/tlr5_synteny.tsv"

for sp in chicken zebrafinch amgo hofi lego pisi; do
  grep -q "^${sp}	" "$OUTDIR/tlr5_synteny.tsv" || die "no synteny rows for $sp"
done

echo
echo "done. Copy to your laptop with:"
echo "  scp -r <user>@<cluster>:$BASE/$OUTDIR <local folder>/"
__END_TLR5_FIGURE_INPUTS_SH__

# ---- tlr5_main_figure.R
cat > "$S/tlr5_main_figure.R" << '__END_TLR5_MAIN_FIGURE_R__'
# =============================================================================
# TLR5 main figure
#   A  reference species tree (cladogram, branch lengths arbitrary) with a
#      single bar marking TLR5 loss of function on the stem of the four finches
#   B  local gene order: one thin line per species with equal-sized solid
#      arrows on it (order only, not to scale), oriented so TLR5 points right;
#      species names at the start of each line, rows level with the tree tips
#   C  the first premature stop shared by all four finches: protein alignment
#      over the corresponding codon-aligned nucleotides
#
# Inputs (working directory): tlr5_synteny.tsv, tlr5_NT.fasta (MACSE NT output)
# Outputs: TLR5_main_figure.pdf / .png
#          tlr5_branch_events.tsv   each distinct disruption and its branch
#          tlr5_species_counts.tsv  disruptions carried by each finch
# Run:   Rscript tlr5_main_figure.R
# Setup: install.packages(c("ggplot2", "patchwork"))
#
# Topology: (chicken, (zebra finch, (house finch, (pine siskin,
#            (American goldfinch, lesser goldfinch))))). A reference tree from
#            the literature (cite it), not inferred here.
# =============================================================================

library(ggplot2)
library(patchwork)

# ---- settings ---------------------------------------------------------------
tip_order <- c("chicken", "zebrafinch", "hofi", "pisi", "amgo", "lego")
tip_label <- c(chicken = "Chicken", zebrafinch = "Zebra finch",
               hofi = "House finch", pisi = "Pine siskin",
               amgo = "American goldfinch", lego = "Lesser goldfinch")
finches <- c("hofi", "pisi", "amgo", "lego")
tip_y   <- setNames(rev(seq_along(tip_order)), tip_order)  # chicken on top
y_lim   <- c(0.4, 6.6)

col_intact <- "#b2182b"     # intact TLR5, loss-of-function bar
col_pseudo <- "#f4a6a6"     # pseudogenized TLR5
col_neigh  <- "gray70"      # neighboring genes
col_stop   <- "#d7301f"     # premature stop / stop-creating base
tree_xmax  <- 13            # right limit of the tree panel; raise it if the
                            # tip names are clipped, lower it to enlarge the tree
col_fs     <- "#4393c3"     # frameshift (drawn as a gap)
flank_codons <- 10          # codons shown either side of the first stop
tag_x_B <- 0.12             # horizontal position of the "B" tag (0 = left edge
                            # of panel B, 1 = right edge); raise to move it right


# ---- alignment and disruptions ----------------------------------------------
read_fasta <- function(path) {
  x  <- readLines(path)
  h  <- grep("^>", x)
  nm <- sub("^>(\\S+).*", "\\1", x[h])
  en <- c(h[-1] - 1, length(x))
  setNames(mapply(function(a, b) toupper(paste(x[(a + 1):b], collapse = "")),
                  h, en), nm)
}
nt <- read_fasta("tlr5_NT.fasta")
stopifnot(all(tip_order %in% names(nt)))
chars <- lapply(tip_order, function(n) strsplit(nt[[n]], "")[[1]])
names(chars) <- tip_order
L <- length(chars[[1]])
stopifnot(all(lengths(chars) == L))

# alignment column -> chicken codon number
ck_codon <- pmax(1, ceiling(cumsum(!chars[["chicken"]] %in% c("-", "!")) / 3))

# frameshifts ('!' runs) and internal stops on MACSE's codon grid
find_events <- function(ch) {
  bang  <- ch == "!"
  fs    <- which(bang & c(TRUE, !head(bang, -1)))
  st    <- seq(1, length(ch) - 2, by = 3)
  cod   <- vapply(st, function(j) paste(ch[j:(j + 2)], collapse = ""), "")
  real  <- which(nchar(gsub("-", "", cod)) > 0)
  stops <- real[cod[real] %in% c("TAA", "TAG", "TGA")]
  stops <- setdiff(stops, tail(real, 1))  # terminal stop is normal
  rbind(data.frame(type = rep("frameshift", length(fs)), col = fs),
        data.frame(type = rep("stop", length(stops)), col = st[stops]))
}
ev <- do.call(rbind, lapply(tip_order, function(n) {
  e <- find_events(chars[[n]]); e$seq <- rep(n, nrow(e)); e
}))

# group events across sequences: stops exact, frameshifts within 3 columns
ev <- ev[order(ev$type, ev$col), ]
first <- numeric(0); ctype <- character(0); ev$cluster <- NA_integer_
for (i in seq_len(nrow(ev))) {
  w   <- if (ev$type[i] == "stop") 0 else 3
  hit <- which(ctype == ev$type[i] & ev$col[i] - first <= w)
  if (length(hit)) {
    ev$cluster[i] <- hit[1]
  } else {
    first <- c(first, ev$col[i]); ctype <- c(ctype, ev$type[i])
    ev$cluster[i] <- length(first)
  }
}

# parsimony placement on the fixed topology: each event goes on the branch
# leading to the smallest clade containing every sequence that carries it
clades <- list(finch_stem = c("hofi", "pisi", "amgo", "lego"),
               spinus     = c("pisi", "amgo", "lego"),
               goldfinch  = c("amgo", "lego"))
assign_branch <- function(s) {
  s <- unique(s)
  if (any(!s %in% finches)) return("in_reference")
  for (b in names(clades)) if (setequal(s, clades[[b]])) return(b)
  if (length(s) == 1) return(s)
  "conflict"
}
br <- do.call(rbind, lapply(split(ev, ev$cluster), function(d)
  data.frame(type = d$type[1], aln_col = min(d$col),
             chicken_codon = ck_codon[min(d$col)],
             carriers = paste(sort(unique(d$seq)), collapse = ","),
             branch = assign_branch(d$seq))))
br <- br[order(br$aln_col), ]

carried <- do.call(rbind, lapply(finches, function(n) {
  e <- ev[ev$seq == n, ]
  data.frame(species = n, stops = sum(e$type == "stop"),
             frameshifts = sum(e$type == "frameshift"))
}))
carried$total <- carried$stops + carried$frameshifts

cat("distinct disruptions per branch:\n")
print(table(branch = br$branch, type = br$type))
cat(sprintf("total distinct: %d stops, %d frameshifts\n",
            sum(br$type == "stop"), sum(br$type == "frameshift")))
cat("disruptions carried by each finch:\n")
print(carried, row.names = FALSE)
if (any(br$branch %in% c("in_reference", "conflict")))
  warning("some events fall in a reference or conflict with the topology")

write.table(br, "tlr5_branch_events.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)
write.table(carried, "tlr5_species_counts.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)


# =============================================================================
# Panel A: tree
# =============================================================================
tip_x <- 7.5
y4 <- mean(tip_y[c("amgo", "lego")])
y3 <- mean(c(tip_y[["pisi"]], y4))
y2 <- mean(c(tip_y[["hofi"]], y3))
y1 <- mean(c(tip_y[["zebrafinch"]], y2))
y0 <- mean(c(tip_y[["chicken"]], y1))

pos <- rbind(
  data.frame(id = c("R", "N1", "N2", "N3", "N4"),
             x  = c(0, 1.5, 3, 4.5, 6), y = c(y0, y1, y2, y3, y4)),
  data.frame(id = tip_order, x = tip_x, y = unname(tip_y[tip_order])))
rownames(pos) <- pos$id

edges <- data.frame(
  parent = c("R", "R", "N1", "N1", "N2", "N2", "N3", "N3", "N4", "N4"),
  child  = c("chicken", "N1", "zebrafinch", "N2", "hofi", "N3",
             "pisi", "N4", "amgo", "lego"))
edges$x0 <- pos[edges$parent, "x"]
edges$x1 <- pos[edges$child, "x"]
edges$y  <- pos[edges$child, "y"]
vert <- do.call(rbind, lapply(split(edges, edges$parent), function(d)
  data.frame(x = d$x0[1], ymin = min(d$y), ymax = max(d$y))))

# single bar on the stem leading to the four finches
if (!any(br$branch == "finch_stem"))
  warning("no disruptions are shared by all four finches")
stem_mid <- (pos["N1", "x"] + pos["N2", "x"]) / 2
loss_bar <- data.frame(xmin = stem_mid - 0.08, xmax = stem_mid + 0.08,
                       ymin = y2 - 0.22, ymax = y2 + 0.22)

tips <- pos[tip_order, ]
tips$label <- tip_label[tips$id]

pTree <- ggplot() +
  geom_segment(aes(x = -0.4, xend = 0, y = y0, yend = y0), linewidth = 0.5) +
  geom_segment(data = edges, aes(x = x0, xend = x1, y = y, yend = y),
               linewidth = 0.5) +
  geom_segment(data = vert, aes(x = x, xend = x, y = ymin, yend = ymax),
               linewidth = 0.5) +
  geom_rect(data = loss_bar,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = col_intact) +
  annotate("text", x = stem_mid - 0.1, y = y2 - 0.3,
           label = "TLR5 loss\nof function", vjust = 1, size = 3.5,
           lineheight = 0.9, color = col_intact, fontface = "bold") +
  geom_text(data = tips, aes(x = x + 0.15, y = y, label = label),
            hjust = 0, size = 3.7, fontface = "bold") +
  scale_x_continuous(limits = c(-0.5, tree_xmax), expand = c(0, 0)) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  theme_void()


# =============================================================================
# Panel B: gene order
# =============================================================================
syn <- read.delim("tlr5_synteny.tsv", stringsAsFactors = FALSE)
stopifnot(all(tip_order %in% syn$species))
syn <- syn[syn$species %in% tip_order, ]

layout_genes <- function(d) {
  d   <- d[order(d$start), ]
  f   <- which(d$focal == 1)[1]
  rel <- seq_len(nrow(d)) - f
  if (d$strand[f] == "-") {               # mirror so TLR5 points right
    rel <- -rel
    d$strand <- ifelse(d$strand == "+", "-", "+")
  }
  d$pos <- rel
  d
}
syn <- do.call(rbind, lapply(split(syn, syn$species), layout_genes))
syn$y <- unname(tip_y[syn$species])
syn$role <- ifelse(syn$focal == 1,
                   ifelse(syn$species %in% finches,
                          "TLR5 (pseudogene)", "TLR5 (intact)"),
                   "neighboring gene")

arrow_poly <- function(p, y, strand, id, half = 0.4, body = 0.13,
                       head = 0.24, headlen = 0.16) {
  x0 <- p - half; x1 <- p + half
  xs <- if (strand == "+") {
    c(x0, x1 - headlen, x1 - headlen, x1, x1 - headlen, x1 - headlen, x0)
  } else {
    c(x1, x0 + headlen, x0 + headlen, x0, x0 + headlen, x0 + headlen, x1)
  }
  ys <- c(y - body, y - body, y - head, y, y + head, y + body, y + body)
  data.frame(id = id, x = xs, y = ys)
}
polys <- do.call(rbind, lapply(seq_len(nrow(syn)), function(i) {
  p <- arrow_poly(syn$pos[i], syn$y[i], syn$strand[i], i)
  p$role <- syn$role[i]
  p
}))

# one line per species, spanning its genes, with the name at the start
lines <- do.call(rbind, lapply(split(syn, syn$species), function(d)
  data.frame(species = d$species[1], y = d$y[1],
             x0 = min(d$pos) - 0.55, x1 = max(d$pos) + 0.55)))
lines$label <- tip_label[lines$species]

role_cols <- c(`TLR5 (intact)` = col_intact, `TLR5 (pseudogene)` = col_pseudo,
               `neighboring gene` = col_neigh)

pSyn <- ggplot() +
  geom_segment(data = lines, aes(x = x0, xend = x1, y = y, yend = y),
               color = "gray40", linewidth = 0.5) +
  geom_polygon(data = polys,
               aes(x, y, group = id, fill = role, color = role),
               linewidth = 0.3) +
  geom_text(data = syn, aes(x = pos, y = y + 0.3, label = gene),
            vjust = 0, size = 3.1) +
  # species name immediately left of the start of its own line
  geom_text(data = lines, aes(x = x0 - 0.15, y = y, label = label),
            hjust = 1, size = 3.7, fontface = "bold") +
  scale_fill_manual(values = role_cols, name = NULL) +
  scale_color_manual(values = role_cols, guide = "none") +
  scale_x_continuous(limits = c(-5.2, 2.7), expand = c(0, 0)) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 10),
        plot.tag.position = c(tag_x_B, 0.97))


# =============================================================================
# Panel C: protein over nucleotides at the first shared stop
# =============================================================================
b <- c("T", "C", "A", "G")
codon_table <- setNames(
  strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
           "")[[1]],
  paste0(rep(b, each = 16), rep(rep(b, each = 4), 4), rep(b, 16)))
translate <- function(cod) {
  if (cod == "---") return("-")
  if (grepl("!", cod, fixed = TRUE)) return("!")   # frameshift -> gap below
  if (grepl("[^ACGT]", cod)) return("X")
  codon_table[[cod]]
}

stop_col <- min(br$aln_col[br$branch == "finch_stem" & br$type == "stop"])
from <- max(1, stop_col - 3 * flank_codons)
to   <- min(L, stop_col + 2 + 3 * flank_codons)
cs   <- seq(from, to - 2, by = 3)           # codon starts in the window

# protein block, BoxShade shading within the window
aa_grp <- c(I = 1, L = 1, M = 1, V = 1, F = 2, W = 2, Y = 2, K = 3, R = 3,
            H = 3, D = 4, E = 4, N = 5, Q = 5, S = 6, T = 6, A = 7, G = 7,
            C = 8, P = 9)
thr <- ceiling(length(tip_order) / 2)

aam <- t(vapply(tip_order, function(n)
  vapply(cs, function(j) translate(paste(chars[[n]][j:(j + 2)], collapse = "")),
         ""), character(length(cs))))
shade_p <- matrix("none", nrow(aam), ncol(aam))
for (j in seq_len(ncol(aam))) {
  colj <- aam[, j]
  res  <- colj[colj %in% names(aa_grp)]
  if (!length(res)) next
  tab  <- sort(table(res), decreasing = TRUE); cons <- names(tab)[1]
  if (tab[1] >= thr) {
    shade_p[colj == cons, j] <- "identical"
    sim <- colj %in% names(aa_grp) & colj != cons & aa_grp[colj] == aa_grp[cons]
    shade_p[which(sim), j] <- "similar"
  }
}
shade_p[aam == "*"] <- "stop"
shade_p[aam == "!"] <- "frameshift"
aam[aam == "!"] <- "-"

prot <- data.frame(seq = rep(tip_order, times = length(cs)),
                   x   = rep(cs + 1, each = length(tip_order)),
                   aa  = as.vector(aam), shade = as.vector(shade_p))
prot$y <- unname(tip_y[prot$seq])

# nucleotide block: consensus shading, stop codon highlighted in carriers,
# with the base(s) differing from the intact zebra finch codon marked
ntm <- t(vapply(tip_order, function(n) chars[[n]][from:to],
                character(to - from + 1)))
shade_n <- matrix("none", nrow(ntm), ncol(ntm))
for (j in seq_len(ncol(ntm))) {
  colj <- ntm[, j]
  res  <- colj[colj %in% c("A", "C", "G", "T")]
  if (!length(res)) next
  tab  <- sort(table(res), decreasing = TRUE)
  if (tab[1] >= thr) shade_n[colj == names(tab)[1], j] <- "identical"
}
shade_n[ntm == "!"] <- "frameshift"
stop_cols_win <- (stop_col:(stop_col + 2)) - from + 1
carriers <- unique(ev$seq[ev$type == "stop" & ev$col == stop_col])
for (n in carriers) {                       # whole stop codon, solid red
  shade_n[match(n, tip_order), stop_cols_win] <- "stop"
}
ntd <- data.frame(seq  = rep(tip_order, times = ncol(ntm)),
                  x    = rep(from:to, each = length(tip_order)),
                  base = as.vector(ntm), shade = as.vector(shade_n))
ntd$y <- unname(tip_y[ntd$seq])

fills <- c(identical = "black", similar = "gray65", none = "white",
           stop = col_stop, frameshift = col_fs)
txt <- c(identical = "white", similar = "black", none = "black",
         stop = "white", frameshift = "white")
fill_breaks <- "stop"                       # the only legend entry
fill_labels <- "Stop codon"

stop_box <- data.frame(xmin = stop_col - 0.5, xmax = stop_col + 2.5,
                       ymin = 0.5, ymax = 6.5)
x_lim <- c(from - 0.5, to + 0.5)
y_axis <- scale_y_continuous(breaks = unname(tip_y),
                             labels = unname(tip_label[names(tip_y)]),
                             limits = c(0.5, 6.5), expand = c(0, 0))
aln_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        plot.title = element_text(size = 11, face = "bold"))

pProt <- ggplot(prot, aes(x, y)) +
  geom_tile(aes(fill = shade), width = 2.9, height = 0.9) +
  geom_text(aes(label = aa, color = shade), size = 3.6, family = "mono",
            fontface = "bold") +
  geom_rect(data = stop_box, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = col_stop, linewidth = 0.8) +
  scale_fill_manual(values = fills, guide = "none") +
  scale_color_manual(values = txt, guide = "none") +
  scale_x_continuous(limits = x_lim, expand = c(0, 0),
                     breaks = cs + 1, labels = ck_codon[cs], position = "top") +
  y_axis +
  labs(title = "Protein", x = NULL, y = NULL) +
  aln_theme

pNuc <- ggplot(ntd, aes(x, y)) +
  geom_tile(aes(fill = shade), width = 0.95, height = 0.9) +
  geom_text(aes(label = base, color = shade), size = 3.1, family = "mono",
            fontface = "bold") +
  geom_rect(data = stop_box, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = col_stop, linewidth = 0.8) +
  scale_fill_manual(values = fills, limits = names(fills),
                    breaks = fill_breaks, labels = fill_labels, name = NULL) +
  scale_color_manual(values = txt, guide = "none") +
  scale_x_continuous(limits = x_lim, expand = c(0, 0), breaks = NULL) +
  y_axis +
  labs(title = "Nucleotide (codon-aligned)", x = NULL, y = NULL) +
  aln_theme +
  theme(legend.position = "bottom", legend.text = element_text(size = 10))

# console summary for the text
for (n in tip_order) {
  cod <- paste(chars[[n]][stop_col:(stop_col + 2)], collapse = "")
  cat(sprintf("%-11s codon %s -> %s\n", n, cod, translate(cod)))
}
cat(sprintf("first shared stop at chicken codon %d\n", ck_codon[stop_col]))


# =============================================================================
# Assemble and save
# =============================================================================
fig <- ((pTree | pSyn) + plot_layout(widths = c(1.5, 2.6))) /
  wrap_elements(full = pProt / pNuc) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A")

pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
ggsave("TLR5_main_figure.pdf", fig, width = 14, height = 10.5, device = pdf_dev, bg = "white")
ggsave("TLR5_main_figure.png", fig, width = 14, height = 10.5, dpi = 300, bg = "white")
cat("wrote TLR5_main_figure.pdf and TLR5_main_figure.png\n")
__END_TLR5_MAIN_FIGURE_R__

# ---- tlr5_figure.R
cat > "$S/tlr5_figure.R" << '__END_TLR5_FIGURE_R__'
# =============================================================================
# TLR5 figure
#   A  local gene order: TLR5 and two neighbors on each side, equal-sized
#      arrows (order only, not to scale), oriented so TLR5 points right;
#      finch TLR5 outlined dashed (pseudogene)
#   B  MACSE protein alignment in black-and-white BoxShade style, wrapped in
#      blocks; frameshift codons shown as gaps; premature stops (red) and
#      frameshift gaps (blue) are the only colored cells
#
# Inputs (in the working directory): tlr5_synteny.tsv, tlr5_AA.fasta
# Run:   Rscript tlr5_figure.R
# Setup: install.packages(c("ggplot2", "gggenes", "patchwork"))
# =============================================================================

library(ggplot2)
library(gggenes)
library(patchwork)

seq_order <- c("chicken", "zebrafinch", "hofi", "pisi", "amgo", "lego")
seq_label <- c(chicken = "Chicken", zebrafinch = "Zebra finch",
               hofi = "House finch", amgo = "American goldfinch",
               lego = "Lesser goldfinch", pisi = "Pine siskin")
finches <- c("hofi", "amgo", "lego", "pisi")
lev <- rev(seq_order)                     # top-to-bottom order on the y axis
lab <- rev(unname(seq_label[seq_order]))

col_tlr5 <- "#b2182b"
col_stop <- "#d7301f"
col_fs   <- "#4393c3"
block_width <- 80                         # residues per alignment block


# -----------------------------------------------------------------------------
# Panel A: gene order
# -----------------------------------------------------------------------------
syn <- read.delim("tlr5_synteny.tsv", stringsAsFactors = FALSE)
stopifnot(all(seq_order %in% syn$species))

layout_genes <- function(d) {
  d   <- d[order(d$start), ]
  f   <- which(d$focal == 1)[1]
  rel <- seq_len(nrow(d)) - f
  if (d$strand[f] == "-") {               # mirror so TLR5 points right
    rel <- -rel
    d$strand <- ifelse(d$strand == "+", "-", "+")
  }
  d$pos  <- rel
  d$xmin <- ifelse(d$strand == "+", rel - 0.4, rel + 0.4)
  d$xmax <- ifelse(d$strand == "+", rel + 0.4, rel - 0.4)
  d
}
syn <- do.call(rbind, lapply(split(syn, syn$species), layout_genes))

syn$role   <- ifelse(syn$focal == 1, "TLR5", "neighbor")
syn$status <- factor(ifelse(syn$focal == 1 & syn$species %in% finches,
                            "pseudogene", "intact"),
                     levels = c("intact", "pseudogene"))
syn$y <- factor(syn$species, levels = lev, labels = lab)

pA <- ggplot(syn, aes(xmin = xmin, xmax = xmax, y = y)) +
  geom_gene_arrow(aes(fill = role, linetype = status), color = "black",
                  arrowhead_height = unit(4, "mm"),
                  arrowhead_width  = unit(2, "mm"),
                  arrow_body_height = unit(3, "mm")) +
  geom_text(aes(x = pos, label = gene), vjust = -1.7, size = 2.4) +
  scale_fill_manual(values = c(TLR5 = col_tlr5, neighbor = "gray80"),
                    guide = "none") +
  scale_linetype_manual(values = c(intact = "solid", pseudogene = "dashed"),
                        name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(-2.6, 2.6)) +
  labs(x = NULL, y = NULL) +
  theme_genes() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())


# -----------------------------------------------------------------------------
# Panel B: protein alignment
# -----------------------------------------------------------------------------
read_fasta <- function(path) {
  x  <- readLines(path)
  h  <- grep("^>", x)
  nm <- sub("^>(\\S+).*", "\\1", x[h])
  en <- c(h[-1] - 1, length(x))
  setNames(mapply(function(a, b) toupper(paste(x[(a + 1):b], collapse = "")),
                  h, en), nm)
}
aa <- read_fasta("tlr5_AA.fasta")
stopifnot(all(seq_order %in% names(aa)))
stopifnot(length(unique(nchar(aa[seq_order]))) == 1)

mat <- do.call(rbind, lapply(seq_order, function(n) strsplit(aa[[n]], "")[[1]]))
rownames(mat) <- seq_order

# MACSE writes frameshifted codons as '!'; draw them as gaps but remember them
fs <- mat == "!"
mat[fs] <- "-"

# internal stops: a '*' followed later by any residue in the same sequence
stop_int <- matrix(FALSE, nrow(mat), ncol(mat))
for (i in seq_len(nrow(mat))) {
  last <- max(which(mat[i, ] != "-"))
  st   <- which(mat[i, ] == "*")
  stop_int[i, st[st < last]] <- TRUE
}

cat("per-sequence check (internal stops / frameshifts):\n")
print(data.frame(row.names = seq_order,
                 internal_stops = rowSums(stop_int),
                 frameshifts    = rowSums(fs)))

# BoxShade-style shading: identical to the column consensus -> black,
# same physicochemical group as the consensus -> gray
aa_grp <- c(I = 1, L = 1, M = 1, V = 1, F = 2, W = 2, Y = 2, K = 3, R = 3,
            H = 3, D = 4, E = 4, N = 5, Q = 5, S = 6, T = 6, A = 7, G = 7,
            C = 8, P = 9)
thr   <- ceiling(nrow(mat) / 2)           # consensus needs >= 3 of 6
shade <- matrix("none", nrow(mat), ncol(mat))
for (j in seq_len(ncol(mat))) {
  col <- mat[, j]
  res <- col[col %in% names(aa_grp)]
  if (!length(res)) next
  tab  <- sort(table(res), decreasing = TRUE)
  cons <- names(tab)[1]
  if (tab[1] >= thr) {
    shade[col == cons, j] <- "identical"
    sim <- col %in% names(aa_grp) & col != cons & aa_grp[col] == aa_grp[cons]
    shade[which(sim), j] <- "similar"
  }
}
shade[stop_int] <- "stop"
shade[fs]       <- "frameshift"

df <- data.frame(seq   = rep(seq_order, times = ncol(mat)),
                 col   = rep(seq_len(ncol(mat)), each = nrow(mat)),
                 aa    = as.vector(mat),
                 shade = as.vector(shade))

# wrap into blocks, labeled by chicken residue numbers
df$block <- (df$col - 1) %/% block_width
df$x     <- (df$col - 1) %% block_width + 1
df$y     <- factor(df$seq, levels = lev, labels = lab)

ck_num <- cumsum(mat["chicken", ] != "-")
blocks <- split(seq_len(ncol(mat)), (seq_len(ncol(mat)) - 1) %/% block_width)
strips <- vapply(seq_along(blocks), function(k) {
  cc <- blocks[[k]]
  sprintf("chicken %d-%d", max(1, ck_num[cc[1]]), ck_num[tail(cc, 1)])
}, "")
df$block_lab <- factor(strips[df$block + 1], levels = unique(strips))
n_blocks <- length(blocks)

fills  <- c(identical = "black", similar = "gray65", none = "white",
            stop = col_stop, frameshift = col_fs)
txtcol <- c(identical = "white", similar = "black", none = "black",
            stop = "white", frameshift = "white")

pB <- ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = shade), width = 0.95, height = 0.9) +
  geom_text(aes(label = aa, color = shade), size = 1.9, family = "mono") +
  scale_fill_manual(values = fills,
                    breaks = c("identical", "similar", "stop", "frameshift"),
                    labels = c("identical (>= 3 of 6)", "similar",
                               "premature stop", "frameshift (gap)"),
                    name = NULL) +
  scale_color_manual(values = txtcol, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~ block_lab, ncol = 1) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(hjust = 0, size = 7),
        legend.position = "top",
        panel.spacing.y = unit(2, "mm"))


# -----------------------------------------------------------------------------
# Assemble and save
# -----------------------------------------------------------------------------
fig_h <- 2.4 + n_blocks * 0.9
fig <- pA / pB +
  plot_layout(heights = c(2.2, n_blocks * 0.9)) +
  plot_annotation(tag_levels = "A")

pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
ggsave("TLR5_figure.pdf", fig, width = 11, height = fig_h, device = pdf_dev, bg = "white")
ggsave("TLR5_figure.png", fig, width = 11, height = fig_h, dpi = 300, bg = "white")
cat("wrote TLR5_figure.pdf and TLR5_figure.png\n")
__END_TLR5_FIGURE_R__

# ---- tlr1_cluster_inputs.sh
cat > "$S/tlr1_cluster_inputs.sh" << '__END_TLR1_CLUSTER_INPUTS_SH__'
#!/bin/bash
# =============================================================================
# Inputs for the TLR1-cluster synteny figure. Run on the cluster from finches/.
#
#   human    T2T-CHM13v2.0 RefSeq annotation (GCF_009914755.1): TLR10, TLR1, TLR6
#            uses the copy in immune_annot/ if it is there, otherwise downloads
#            it with the NCBI `datasets` tool
#   chicken  already on disk: gaga_genomic.gff and gaga_prot.faa
#            (GCF_016700215.2, paternal white leghorn; symbols TLR1A / TLR1B)
#   birds    the miniprot GFFs (chicken proteins aligned to each genome) of the
#            finch primaries and the zebra finch assembly, all handled the same
#            way: annotation/miniprot/*/*_pri.gff and bTaeGut7_mat.gff
#
# Then runs scripts/tlr1_cluster_synteny.py, which writes
#   tlr1_synteny/tlr1_synteny.tsv
#
# Override locations with environment variables if needed, e.g.
#   HUMAN_LOCAL=/path/to/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz bash ...
# =============================================================================

set -euo pipefail

BASE=${BASE:?set BASE to the project directory}
SCRIPTS=${SCRIPTS:-$BASE/scripts}
HUMAN_LOCAL=${HUMAN_LOCAL:-}   # optional local copy of GCF_009914755.1 genomic.gff.gz; downloaded if empty
HUMAN_ACC=GCF_009914755.1

die() { echo "ERROR: $*" >&2; exit 1; }

cd "$BASE" || die "cannot cd to $BASE"
mkdir -p refs/human tlr1_synteny

# ---- chicken: must already be here ------------------------------------------
for f in gaga_genomic.gff gaga_prot.faa; do
  [ -s "$f" ] || die "missing $f (the paternal chicken annotation and proteins)"
done
echo "chicken: gaga_genomic.gff ($(awk -F'\t' '$3=="gene"' gaga_genomic.gff | wc -l) genes)"

# ---- human: local copy first, download otherwise ----------------------------
if [ -s refs/human/genomic.gff.gz ]; then
  echo "human: using refs/human/genomic.gff.gz"
elif [ -s "$HUMAN_LOCAL" ]; then
  cp "$HUMAN_LOCAL" refs/human/genomic.gff.gz
  echo "human: copied $HUMAN_LOCAL"
else
  command -v datasets >/dev/null 2>&1 \
    || die "no local human GFF and the NCBI 'datasets' tool is not on PATH (activate the environment you used for the zebra finch download)"
  tmp=$(mktemp -d)
  ( cd "$tmp" && datasets download genome accession "$HUMAN_ACC" \
      --include gff3,protein --filename human.zip && unzip -q human.zip )
  d="$tmp/ncbi_dataset/data/$HUMAN_ACC"
  [ -s "$d/genomic.gff" ] || die "download did not produce genomic.gff"
  gzip -c "$d/genomic.gff" > refs/human/genomic.gff.gz
  [ -s "$d/protein.faa" ] && cp "$d/protein.faa" refs/human/protein.faa
  rm -rf "$tmp"
  echo "human: downloaded $HUMAN_ACC"
fi

# the three genes that define the human cluster must be annotated
for g in TLR10 TLR1 TLR6; do
  zcat refs/human/genomic.gff.gz | awk -F'\t' -v g="gene=$g;" '$3=="gene" && index($9,g)' \
    | head -1 | cut -f1,4,5,7 | sed "s/^/human $g: /" \
    | grep . || die "human $g not found in the annotation"
done

# optional: human proteins (not needed for the figure, useful for later trees)
[ -s refs/human/protein.faa ] || echo "note: refs/human/protein.faa not present (only needed for trees/RBH)"

# ---- build the table ---------------------------------------------------------
[ -s "$SCRIPTS/tlr1_cluster_synteny.py" ] || die "missing $SCRIPTS/tlr1_cluster_synteny.py"
python3 "$SCRIPTS/tlr1_cluster_synteny.py" --out tlr1_synteny/tlr1_synteny.tsv

echo
echo "next, on your laptop:"
echo "  scp <user>@<cluster>:$BASE/tlr1_synteny/tlr1_synteny.tsv <local folder>/tlr5/"
echo "  cd <local folder>/tlr5 && Rscript tlr1_synteny_figure.R"
__END_TLR1_CLUSTER_INPUTS_SH__

# ---- tlr1_cluster_synteny.py
cat > "$S/tlr1_cluster_synteny.py" << '__END_TLR1_CLUSTER_SYNTENY_PY__'
#!/usr/bin/env python3
"""
Gene order around the TLR1-family cluster in human, chicken, zebra finch and
the finches.

Cluster genes (the genes that define the locus in each species):
  human (CHM13)  TLR10, TLR1, TLR6
  chicken        TLR1A, TLR1B           (NCBI symbols)
  zebra finch    the miniprot loci named TLR1A and TLR1B (same as the finches)
  finches        the miniprot loci named TLR1A and TLR1B

For each species the script finds the cluster, adds --flank protein-coding
genes beyond each end, and writes one table row per gene. Every species is
oriented against human using the flanking genes they share (a species is
flipped if the shared genes run in the opposite order), so the tables can be
drawn side by side. Each species keeps its own gene symbols.

Human and chicken use their NCBI annotations, read as gene features with
gene_biotype=protein_coding. The zebra finch and the finches are all treated
the same way, from miniprot alignments of chicken proteins to their genomes:
overlapping hits on the same strand are merged into loci, and each locus is
named after the chicken gene of its best-scoring hit, so their genes carry
chicken symbols (the zebra finch's own annotation is not used).

Usage (from finches/):
  python3 scripts/tlr1_cluster_synteny.py
  python3 scripts/tlr1_cluster_synteny.py --flank 4 --out tlr1_synteny/tlr1_synteny.tsv

Output columns: species, order, gene, contig, start, end, strand, cluster
  order   1..n, left to right after orientation
  strand  after orientation (flipped species have their strands swapped)
  cluster 1 for the cluster genes
Only the Python standard library is needed.
"""

import argparse
import glob
import gzip
import os
import re
import sys
from collections import defaultdict
from itertools import combinations

HUMAN_CLUSTER = {"TLR10", "TLR1", "TLR6"}
CHICKEN_CLUSTER = {"TLR1A", "TLR1B"}
FINCH_ORDER = ["hofi", "pisi", "amgo", "lego"]


def die(msg):
    sys.exit(f"ERROR: {msg}")


def attrs(s):
    return dict(kv.split("=", 1) for kv in s.split(";") if "=" in kv)


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


# ----------------------------------------------------------------------------
# readers
# ----------------------------------------------------------------------------
def read_ncbi_genes(path):
    genes = []
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            a = attrs(c[8])
            if a.get("gene_biotype") != "protein_coding":
                continue
            sym = a.get("gene") or a.get("Name")
            if sym:
                genes.append({"chr": c[0], "s": int(c[3]), "e": int(c[4]),
                              "strand": c[6], "gene": sym})
    return genes


def read_chicken_map(path):
    """protein_id -> gene symbol from the chicken CDS lines."""
    p2g = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "CDS":
                continue
            a = attrs(c[8])
            if "protein_id" in a and "gene" in a:
                p2g[a["protein_id"]] = a["gene"]
    return p2g


def read_lengths(path):
    lens, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            lens[name] = 0
        elif name:
            lens[name] += len(line.strip())
    return lens


def read_miniprot_loci(path, p2g, plen, min_pos, min_cov):
    hits = []
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "mRNA":
            continue
        a = attrs(c[8])
        if float(a.get("Positive", 0)) < min_pos:
            continue
        tgt = a["Target"].split()
        prot, qs, qe = tgt[0], int(tgt[1]), int(tgt[2])
        if prot in plen and (qe - qs + 1) / plen[prot] < min_cov:
            continue
        hits.append({"chr": c[0], "s": int(c[3]), "e": int(c[4]),
                     "strand": c[6],
                     "score": float(c[5]) if c[5] != "." else 0.0,
                     "gene": p2g.get(prot, prot)})
    hits.sort(key=lambda h: (h["chr"], h["strand"], h["s"]))
    loci, cur = [], None
    for h in hits:
        if (cur and h["chr"] == cur["chr"] and h["strand"] == cur["strand"]
                and h["s"] <= cur["e"]):
            cur["e"] = max(cur["e"], h["e"])
            if h["score"] > cur["score"]:
                cur["score"], cur["gene"] = h["score"], h["gene"]
        else:
            if cur:
                loci.append(cur)
            cur = dict(h)
    if cur:
        loci.append(cur)
    return loci


# ----------------------------------------------------------------------------
# neighborhood and orientation
# ----------------------------------------------------------------------------
def neighborhood(genes, cluster, flank):
    """Cluster span plus `flank` genes beyond each end, in genomic order."""
    found = [g for g in genes if g["gene"] in cluster]
    if not found:
        return None, f"none of {sorted(cluster)} found"
    by_chr = defaultdict(list)
    for g in found:
        by_chr[g["chr"]].append(g)
    chrom = max(by_chr, key=lambda k: len(by_chr[k]))
    msgs = []
    if len(by_chr) > 1:
        msgs.append("cluster genes on several contigs: "
                    + ", ".join(f"{k}:{len(v)}" for k, v in by_chr.items())
                    + f"; using {chrom}")
    missing = cluster - {g["gene"] for g in by_chr[chrom]}
    if missing:
        msgs.append(f"missing on {chrom}: {', '.join(sorted(missing))}")
    span_s = min(g["s"] for g in by_chr[chrom])
    span_e = max(g["e"] for g in by_chr[chrom])
    on = sorted((g for g in genes if g["chr"] == chrom),
                key=lambda g: (g["s"], g["e"]))
    left = [g for g in on if g["e"] < span_s][-flank:]
    inside = [g for g in on if g["e"] >= span_s and g["s"] <= span_e]
    right = [g for g in on if g["s"] > span_e][:flank]
    rows = []
    for g in left + inside + right:
        r = dict(g)
        r["cluster"] = 1 if g["gene"] in cluster else 0
        rows.append(r)
    return rows, "; ".join(msgs)


def flip(rows):
    out = []
    for r in reversed(rows):
        r = dict(r)
        r["strand"] = "-" if r["strand"] == "+" else "+"
        out.append(r)
    return out


def orient_to(rows, ref_rank):
    """Flip rows if the non-cluster genes shared with the reference run in
    the opposite order. Returns (rows, note)."""
    shared = [(i, ref_rank[r["gene"]]) for i, r in enumerate(rows)
              if not r["cluster"] and r["gene"] in ref_rank]
    if len(shared) < 2:
        return rows, f"{len(shared)} shared flanking gene(s): not oriented"
    conc = disc = 0
    for (i1, a), (i2, b) in combinations(shared, 2):
        if (i2 - i1) * (b - a) > 0:
            conc += 1
        else:
            disc += 1
    if disc > conc:
        return flip(rows), f"flipped ({len(shared)} shared flanking genes)"
    return rows, f"kept ({len(shared)} shared flanking genes)"


def describe(rows):
    first = next(i for i, r in enumerate(rows) if r["cluster"])
    last = max(i for i, r in enumerate(rows) if r["cluster"])
    s = [f"{r['gene']}({r['strand']})" for r in rows]
    s.insert(last + 1, "]")
    s.insert(first, "[")
    return " ".join(s)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--human-gff", default="refs/human/genomic.gff.gz")
    ap.add_argument("--chicken-gff", default="gaga_genomic.gff")
    ap.add_argument("--chicken-faa", default="gaga_prot.faa")
    ap.add_argument("--miniprot-dir", default="annotation/miniprot")
    ap.add_argument("--flank", type=int, default=3,
                    help="genes shown beyond each end of the cluster")
    ap.add_argument("--min-positive", type=float, default=0.5)
    ap.add_argument("--min-cov", type=float, default=0.5,
                    help="fraction of the chicken protein a finch hit must cover")
    ap.add_argument("--out", default="tlr1_synteny/tlr1_synteny.tsv")
    args = ap.parse_args()

    for f in (args.human_gff, args.chicken_gff, args.chicken_faa):
        if not os.path.isfile(f):
            die(f"not found: {f}")

    tables = {}          # species -> rows (genomic orientation)
    notes = {}

    human = read_ncbi_genes(args.human_gff)
    rows, msg = neighborhood(human, HUMAN_CLUSTER, args.flank)
    if rows is None:
        die(f"human: {msg}")
    tables["human"], notes["human"] = rows, msg

    chicken = read_ncbi_genes(args.chicken_gff)
    rows, msg = neighborhood(chicken, CHICKEN_CLUSTER, args.flank)
    if rows is None:
        die(f"chicken: {msg}")
    tables["chicken"], notes["chicken"] = rows, msg

    p2g = read_chicken_map(args.chicken_gff)
    plen = read_lengths(args.chicken_faa)
    gffs = {}
    for f in sorted(glob.glob(os.path.join(args.miniprot_dir, "*", "*.gff"))):
        base = os.path.basename(f)[:-4]
        if base.startswith("bTaeGut7"):
            gffs["zebrafinch"] = f              # zebra finch assembly
        elif base.endswith("_alt"):
            continue                            # primaries only
        else:
            gffs[re.sub(r"_pri$", "", base)] = f
    if "zebrafinch" in gffs:
        print(f"zebra finch: using {gffs['zebrafinch']} "
              "(chicken proteins aligned with miniprot, not its own annotation)")
    else:
        print("WARNING: no miniprot GFF for the zebra finch (bTaeGut7*.gff) found")
    for sp in [s for s in ["zebrafinch"] + FINCH_ORDER if s in gffs] + \
              sorted(s for s in gffs if s not in ["zebrafinch"] + FINCH_ORDER):
        loci = read_miniprot_loci(gffs[sp], p2g, plen, args.min_positive,
                                  args.min_cov)
        rows, msg = neighborhood(loci, CHICKEN_CLUSTER, args.flank)
        if rows is None:
            print(f"WARNING {sp}: {msg}")
            continue
        tables[sp], notes[sp] = rows, msg

    # orient everything against human --------------------------------------
    human_flank = [r["gene"] for r in tables["human"] if not r["cluster"]]
    ref_rank = {}
    for i, g in enumerate(human_flank):
        ref_rank.setdefault(g, i)

    final = {"human": tables["human"]}
    order = ["human", "chicken", "zebrafinch"] + \
            [s for s in FINCH_ORDER if s in tables] + \
            sorted(s for s in tables if s not in
                   ["human", "chicken", "zebrafinch"] + FINCH_ORDER)
    for sp in order:
        if sp == "human" or sp not in tables:
            continue
        final[sp], note = orient_to(tables[sp], ref_rank)
        notes[sp] = (notes[sp] + "; " if notes[sp] else "") + note

    # write -----------------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as out:
        out.write("\t".join(["species", "order", "gene", "contig", "start",
                             "end", "strand", "cluster"]) + "\n")
        for sp in order:
            if sp not in final:
                continue
            for i, r in enumerate(final[sp], 1):
                out.write("\t".join([sp, str(i), r["gene"], r["chr"],
                                     str(r["s"]), str(r["e"]), r["strand"],
                                     str(r["cluster"])]) + "\n")

    print("gene order (cluster in [ ]; strands after orientation):\n")
    for sp in order:
        if sp not in final:
            continue
        chrom = final[sp][0]["chr"]
        print(f"{sp:11s} {chrom}")
        print(f"            {describe(final[sp])}")
        if notes[sp]:
            print(f"            note: {notes[sp]}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
__END_TLR1_CLUSTER_SYNTENY_PY__

# ---- tlr1_synteny_figure.R
cat > "$S/tlr1_synteny_figure.R" << '__END_TLR1_SYNTENY_FIGURE_R__'
# =============================================================================
# Gene order around the TLR1-family cluster (human TLR10-TLR1-TLR6, bird
# TLR1A/TLR1B): one thin line per species with equal-sized solid arrows on it,
# arrow direction = strand, species name at the start of its line. Order only,
# not to scale. Rows are aligned so every cluster starts at the same x.
#
# Input (working directory): tlr1_synteny.tsv from tlr1_cluster_synteny.py
# Run:    Rscript tlr1_synteny_figure.R
# Setup:  install.packages("ggplot2")
#
# Coloring (set `pairing` below):
#   "velova"  TLR10 and TLR1A share one color; TLR1, TLR6 and TLR1B share
#             another. This is the pairing proposed by Velova et al. (2018);
#             the figure then shows whether position agrees with it.
#   "none"    every cluster gene gets one color, so nothing is presupposed.
# The zebra finch and the finches are all chicken proteins aligned to the
# genome (miniprot), so their genes carry chicken symbols. Any cluster gene
# with another symbol is shown as unassigned in "velova" mode.
# =============================================================================

library(ggplot2)

infile  <- "tlr1_synteny.tsv"
pairing <- "velova"                       # "velova" or "none"

species_order <- c("human", "chicken", "zebrafinch", "hofi", "pisi", "amgo", "lego")
species_label <- c(human = "Human (CHM13)", chicken = "Chicken",
                   zebrafinch = "Zebra finch", hofi = "House finch",
                   pisi = "Pine siskin", amgo = "American goldfinch",
                   lego = "Lesser goldfinch")

col_flank <- "gray72"
group_cols <- c(`TLR10 / TLR1A`        = "#e08214",
                `TLR1 / TLR6 / TLR1B`  = "#2166ac",
                unassigned             = "#8c6bb1",
                `TLR1 family`          = "#b2182b",
                `flanking gene`        = col_flank)

syn <- read.delim(infile, stringsAsFactors = FALSE)
species_order <- species_order[species_order %in% syn$species]
stopifnot(length(species_order) > 0)
syn <- syn[syn$species %in% species_order, ]

y_of  <- setNames(rev(seq_along(species_order)), species_order)
syn$y <- unname(y_of[syn$species])

# align the start of the cluster at x = 1 in every species
first_cl <- tapply(ifelse(syn$cluster == 1, syn$order, NA), syn$species,
                   min, na.rm = TRUE)
syn$x <- syn$order - unname(first_cl[syn$species]) + 1

group_of <- function(gene, cluster) {
  g <- rep("flanking gene", length(gene))
  inc <- cluster == 1
  if (pairing == "none") {
    g[inc] <- "TLR1 family"
  } else {
    g[inc] <- "unassigned"
    g[inc & gene %in% c("TLR10", "TLR1A")]          <- "TLR10 / TLR1A"
    g[inc & gene %in% c("TLR1", "TLR6", "TLR1B")]   <- "TLR1 / TLR6 / TLR1B"
  }
  g
}
syn$group <- group_of(syn$gene, syn$cluster)

arrow_poly <- function(p, y, strand, id, half = 0.4, body = 0.13,
                       head = 0.24, headlen = 0.16) {
  x0 <- p - half; x1 <- p + half
  xs <- if (strand == "+") {
    c(x0, x1 - headlen, x1 - headlen, x1, x1 - headlen, x1 - headlen, x0)
  } else {
    c(x1, x0 + headlen, x0 + headlen, x0, x0 + headlen, x0 + headlen, x1)
  }
  ys <- c(y - body, y - body, y - head, y, y + head, y + body, y + body)
  data.frame(id = id, x = xs, y = ys)
}
polys <- do.call(rbind, lapply(seq_len(nrow(syn)), function(i) {
  p <- arrow_poly(syn$x[i], syn$y[i], syn$strand[i], i)
  p$group <- syn$group[i]
  p
}))

lines <- do.call(rbind, lapply(split(syn, syn$species), function(d)
  data.frame(species = d$species[1], y = d$y[1],
             x0 = min(d$x) - 0.55, x1 = max(d$x) + 0.55)))
lines$label <- species_label[lines$species]

present <- names(group_cols)[names(group_cols) %in% syn$group]
syn$bold <- ifelse(syn$cluster == 1, "bold", "plain")

p <- ggplot() +
  geom_segment(data = lines, aes(x = x0, xend = x1, y = y, yend = y),
               color = "gray40", linewidth = 0.5) +
  geom_polygon(data = polys,
               aes(x, y, group = id, fill = group, color = group),
               linewidth = 0.3) +
  geom_text(data = syn, aes(x = x, y = y + 0.32, label = gene, fontface = bold),
            vjust = 0, size = 3.1) +
  geom_text(data = lines, aes(x = x0 - 0.15, y = y, label = label),
            hjust = 1, size = 3.7, fontface = "bold") +
  scale_fill_manual(values = group_cols, breaks = present, name = NULL) +
  scale_color_manual(values = group_cols, guide = "none") +
  scale_x_continuous(limits = c(min(syn$x) - 3.2, max(syn$x) + 0.8),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, length(species_order) + 0.8),
                     expand = c(0, 0)) +
  labs(caption = paste("Gene order around the TLR1-family cluster; arrow direction shows strand.",
                       "Species are oriented against human using shared flanking genes.")) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 10),
        plot.caption = element_text(size = 8, hjust = 0, color = "gray30"),
        plot.margin = margin(8, 8, 8, 8))

h <- 1.4 + 0.75 * length(species_order)
ggsave("TLR1_synteny.pdf", p, width = 13, height = h, bg = "white")
ggsave("TLR1_synteny.png", p, width = 13, height = h, dpi = 300, bg = "white")
cat("wrote TLR1_synteny.pdf and TLR1_synteny.png\n")
__END_TLR1_SYNTENY_FIGURE_R__

# ---- tlr1_synteny_scaled.R
cat > "$S/tlr1_synteny_scaled.R" << '__END_TLR1_SYNTENY_SCALED_R__'
# =============================================================================
# TLR1-family cluster drawn to scale.
#
# Each gene is a narrow vertical rectangle placed at its real position, so the
# spacing between genes is real (the arrow figure shows order only). Rows are
# centered so that x = 0 is the middle of the cluster (human TLR1; midway
# between TLR1A and TLR1B in the birds), and species are oriented against
# human, exactly as in tlr1_synteny.tsv.
#
# Input (working directory): tlr1_synteny.tsv from tlr1_cluster_synteny.py
# Run:    Rscript tlr1_synteny_scaled.R
# Setup:  install.packages("ggplot2")
#
# Settings below:
#   pairing        "velova" colors TLR10/TLR1A together and TLR1/TLR6/TLR1B
#                  together; "none" gives every cluster gene one color
#   gene_style     "tick": each gene is a narrow rectangle of fixed width at its
#                  midpoint (gene length is not drawn); "span": each rectangle
#                  covers the gene's true start-to-end extent
#   strand_offset  TRUE draws plus-strand genes above the line and minus-strand
#                  genes below it; FALSE centers every gene on the line
#   min_width_kb   narrowest rectangle drawn in "span" style
#   stop_at_gene   the left edge of the window sits at this gene (default PGM2,
#                  the outermost flanking gene in the birds); genes further
#                  left, e.g. human PTTG2 and TBC1D1, are cropped and listed in
#                  the caption. NA shows everything.
#
# Human and chicken coordinates are NCBI gene spans. The zebra finch and the
# finches are chicken proteins aligned with miniprot, so their rectangles are
# the aligned extent of each gene and approximate its true span.
# =============================================================================

library(ggplot2)

infile        <- "tlr1_synteny.tsv"
pairing       <- "velova"          # "velova" or "none"
gene_style    <- "tick"           # "tick" or "span"
strand_offset <- TRUE
min_width_kb  <- 0.5
stop_at_gene  <- "PGM2"             # NA to show all genes

species_order <- c("human", "chicken", "zebrafinch", "hofi", "pisi", "amgo", "lego")
species_label <- c(human = "Human (CHM13)", chicken = "Chicken",
                   zebrafinch = "Zebra finch", hofi = "House finch",
                   pisi = "Pine siskin", amgo = "American goldfinch",
                   lego = "Lesser goldfinch")

col_flank  <- "gray72"
group_cols <- c(`TLR10 / TLR1A`       = "#e08214",
                `TLR1 / TLR6 / TLR1B` = "#2166ac",
                unassigned            = "#8c6bb1",
                `TLR1 family`         = "#b2182b",
                `flanking gene`       = col_flank)

syn <- read.delim(infile, stringsAsFactors = FALSE)
species_order <- species_order[species_order %in% syn$species]
stopifnot(length(species_order) > 0)
syn <- syn[syn$species %in% species_order, ]

# ---- coordinates, left to right in the oriented order, in kb -----------------
parts <- lapply(split(syn, syn$species), function(d) {
  d <- d[order(d$order), ]
  flipped <- nrow(d) > 1 && d$start[1] > d$start[nrow(d)]
  if (flipped) {                         # species was flipped to face human
    xs <- -d$end;   xe <- -d$start
  } else {
    xs <- d$start;  xe <- d$end
  }
  d$xs <- xs / 1000
  d$xe <- xe / 1000
  d
})
syn <- do.call(rbind, parts)

syn$mid <- (syn$xs + syn$xe) / 2

# x = 0 is the middle of the cluster: the median of the cluster genes' midpoints
# (the middle gene, human TLR1, or midway between the two birds' genes)
center <- tapply(ifelse(syn$cluster == 1, syn$mid, NA), syn$species,
                 median, na.rm = TRUE)
syn$xs  <- syn$xs  - unname(center[syn$species])
syn$xe  <- syn$xe  - unname(center[syn$species])
syn$mid <- syn$mid - unname(center[syn$species])

# stop at PGM2: the window starts there and genes further left are cropped
cropped <- syn[0, ]
if (!is.na(stop_at_gene)) {
  hit <- syn$gene == stop_at_gene
  if (any(hit)) {
    left_edge <- min(syn$xs[hit]) - 10
    outside <- if (gene_style == "tick") syn$mid < left_edge else syn$xe < left_edge
    cropped <- syn[outside, ]
    syn <- syn[!outside, ]
    if (nrow(cropped))
      message("cropped left of ", stop_at_gene, ": ",
              paste(cropped$species, cropped$gene, collapse = ", "))
  } else {
    warning(stop_at_gene, " not found in the table; nothing cropped")
  }
}
crop_txt <- if (nrow(cropped))
  paste0("The window starts at ", stop_at_gene, "; further left and not shown: ",
         paste(sprintf("%s %s", species_label[cropped$species], cropped$gene),
               collapse = ", "), ".\n") else ""

if (gene_style == "tick") {
  # fixed narrow width, scaled to the plotted range so it looks the same
  tick_kb <- 0.006 * (max(syn$xe) - min(syn$xs))
  syn$xs  <- syn$mid - tick_kb / 2
  syn$xe  <- syn$mid + tick_kb / 2
} else {
  short <- (syn$xe - syn$xs) < min_width_kb      # keep very short genes visible
  syn$xs[short] <- syn$mid[short] - min_width_kb / 2
  syn$xe[short] <- syn$mid[short] + min_width_kb / 2
}

lim_kb <- 1.04 * max(abs(c(syn$xs, syn$xe)))     # symmetric about x = 0

# ---- rows ---------------------------------------------------------------------
bh      <- if (gene_style == "tick") 0.34 else 0.16     # rectangle height
row_gap <- if (strand_offset) 2 * (bh + 0.47) + 0.3 else 1.5
y_of  <- setNames((rev(seq_along(species_order)) - 1) * row_gap + 1, species_order)
syn$y <- unname(y_of[syn$species])

syn$side <- if (strand_offset) ifelse(syn$strand == "+", 1, -1) else 1
if (strand_offset) {
  syn$ymin <- ifelse(syn$side == 1, syn$y + 0.04, syn$y - 0.04 - bh)
} else {
  syn$ymin <- syn$y - bh / 2
}
syn$ymax <- syn$ymin + bh

# ---- colors ------------------------------------------------------------------
g <- rep("flanking gene", nrow(syn))
inc <- syn$cluster == 1
if (pairing == "none") {
  g[inc] <- "TLR1 family"
} else {
  g[inc] <- "unassigned"
  g[inc & syn$gene %in% c("TLR10", "TLR1A")]         <- "TLR10 / TLR1A"
  g[inc & syn$gene %in% c("TLR1", "TLR6", "TLR1B")]  <- "TLR1 / TLR6 / TLR1B"
}
syn$group <- g
present <- names(group_cols)[names(group_cols) %in% syn$group]

# ---- labels, staggered on two levels per side so neighbors do not collide ----
syn <- syn[order(syn$species, syn$xs), ]
syn$lvl <- ave(seq_along(syn$gene), syn$species, syn$side,
               FUN = function(i) (seq_along(i) - 1) %% 2)
syn$ly <- ifelse(syn$side == 1,
                 syn$ymax + 0.07 + 0.2 * syn$lvl,
                 syn$ymin - 0.07 - 0.2 * syn$lvl)
syn$vj   <- ifelse(syn$side == 1, 0, 1)
syn$face <- ifelse(syn$cluster == 1, "bold", "plain")

spans <- do.call(rbind, lapply(split(syn, syn$species), function(d)
  data.frame(species = d$species[1], y = d$y[1],
             x0 = min(d$xs), x1 = max(d$xe))))

p <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
  geom_segment(data = spans, aes(x = x0, xend = x1, y = y, yend = y),
               color = "gray45", linewidth = 0.4) +
  geom_rect(data = syn,
            aes(xmin = xs, xmax = xe, ymin = ymin, ymax = ymax,
                fill = group, color = group), linewidth = 0.3) +
  geom_text(data = syn,
            aes(x = mid, y = ly, label = gene, fontface = face, vjust = vj),
            size = 2.9) +
  scale_fill_manual(values = group_cols, breaks = present, name = NULL) +
  scale_color_manual(values = group_cols, guide = "none") +
  scale_y_continuous(breaks = unname(y_of),
                     labels = unname(species_label[names(y_of)]),
                     limits = c(min(y_of) - (bh + 0.9), max(y_of) + (bh + 0.9)),
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(-lim_kb, lim_kb), expand = c(0, 0)) +
  labs(x = "Distance from the middle of the cluster (kb)", y = NULL,
       caption = paste0(
         if (gene_style == "tick")
           "Each gene is marked at its midpoint (gene lengths are not drawn); spacing is to scale.\n"
         else "Genes are drawn at their genomic positions and lengths.\n",
         "x = 0 is the middle of the cluster: human TLR1, and midway between TLR1A and TLR1B in the birds.\n",
         "Species are oriented against human using shared flanking genes",
         if (strand_offset) "; plus-strand genes sit above the line, minus-strand genes below.\n"
         else ".\n",
         crop_txt,
         "Zebra finch and finch positions come from chicken proteins aligned with miniprot, ",
         "so they approximate the true gene positions.")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold", size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        legend.position = "bottom", legend.text = element_text(size = 10),
        plot.caption = element_text(size = 8, hjust = 0, color = "gray30"),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(8, 12, 8, 8))

h <- 1.8 + 0.95 * length(species_order)
ggsave("TLR1_synteny_scaled.pdf", p, width = 13, height = h, bg = "white")
ggsave("TLR1_synteny_scaled.png", p, width = 13, height = h, dpi = 300, bg = "white")
cat("wrote TLR1_synteny_scaled.pdf and TLR1_synteny_scaled.png\n")
__END_TLR1_SYNTENY_SCALED_R__


# =============================================================================
# Section 3: check
# =============================================================================
if want check; then
  log "check: inputs, scripts, environments"
  for f in gaga_prot.faa gaga_genomic.gff "${GENOMES[@]}" \
           "$TAGU/protein.faa" "$TAGU/cds_from_genomic.fna" "$TAGU/genomic.gff"; do
    [ -s "$f" ] || die "missing input: $f"
  done
  for f in tlr_flanks.py tlr_similarity_table.py tlr5_analysis.sh tlr5_macse.sh \
           shared_disruptions.py tlr5_figure_inputs.sh tlr5_main_figure.R tlr5_figure.R \
           tlr1_cluster_inputs.sh tlr1_cluster_synteny.py tlr1_synteny_figure.R tlr1_synteny_scaled.R; do
    [ -s "$S/$f" ] || die "missing script: scripts/$f"
  done
  grep -q '^conda activate miniprot' "$S/tlr5_analysis.sh" \
    || die "scripts/tlr5_analysis.sh lacks 'conda activate miniprot' (use the pipeline copy)"
  envs=$(conda env list | awk '!/^#/ && NF {print $1}')
  for e in miniprot exonerate_env macse; do
    grep -qx "$e" <<< "$envs" || die "conda environment '$e' not found"
  done
  echo "all inputs, scripts and environments found"
fi


# =============================================================================
# Section 4: refs
# =============================================================================
if want refs; then
  log "refs: chicken CDS into lab storage"
  mkdir -p refs/chicken
  if [ ! -s refs/chicken/cds_from_genomic.fna ]; then
    [ -s "$CHICKEN_CDS_SRC" ] || die "chicken CDS not found in scratch: $CHICKEN_CDS_SRC"
    cp "$CHICKEN_CDS_SRC" refs/chicken/cds_from_genomic.fna
  fi
  echo "refs/chicken/cds_from_genomic.fna $(md5sum < refs/chicken/cds_from_genomic.fna | cut -d' ' -f1)"
fi


# =============================================================================
# Section 5: proteome
# =============================================================================
if want proteome; then
  log "proteome: chicken proteins mapped genome-wide with miniprot"
  for g in "${GENOMES[@]}"; do
    name=$(basename "$g" .fna); sp=${name%%_*}
    out=annotation/miniprot/$sp/$name.gff
    if [ -s "$out" ]; then
      echo "$name: reusing $out"
      continue
    fi
    mkdir -p "annotation/miniprot/$sp"
    conda activate miniprot
    echo "$name: running miniprot ($(miniprot --version))"
    miniprot -t "$THREADS" -I --gff "$g" gaga_prot.faa > "$out.tmp" 2> "annotation/miniprot/$sp/$name.log"
    mv "$out.tmp" "$out"
    conda deactivate
  done
fi


# =============================================================================
# Section 6: flanks
# =============================================================================
if want flanks; then
  log "flanks: TLR loci and neighbors in every genome"
  mkdir -p annotation/tlr_synteny
  for g in annotation/miniprot/*/*.gff; do
    n=$(basename "$g" .gff)
    python3 "$S/tlr_flanks.py" "$g" gaga_genomic.gff "$n" > "annotation/tlr_synteny/$n.tlr_flanks.tsv"
  done
  cat annotation/tlr_synteny/*.tlr_flanks.tsv | awk 'NR==1 || $1!="sample"' \
    > annotation/tlr_synteny/all_tlr_flanks.tsv
  awk '$3==0{print $1, $8}' annotation/tlr_synteny/all_tlr_flanks.tsv | sort | uniq -c
fi


# =============================================================================
# Section 7: complement
# =============================================================================
if want complement; then
  log "complement: TLR identity/similarity/coverage/copy table"
  python3 "$S/tlr_similarity_table.py" --pri-only --outdir tlr_table
fi


# =============================================================================
# Section 8: tlr5
# =============================================================================
if want tlr5; then
  log "tlr5: locate, translate against chicken, shared stops, exonerate"
  bash "$S/tlr5_analysis.sh"
  [ -s tlr5/results/tlr5_stops.txt ] || die "tlr5 step produced no tlr5/results/tlr5_stops.txt"
fi


# =============================================================================
# Section 9: codon
# =============================================================================
if want codon; then
  log "codon: MACSE codon alignment and shared disruptions"
  [ -s refs/chicken/cds_from_genomic.fna ] || die "run the refs step first"
  bash "$S/tlr5_macse.sh"
  [ -s macse/out/tlr5_NT.fasta ] || die "codon step produced no macse/out/tlr5_NT.fasta"
fi


# =============================================================================
# Section 10: figures
# =============================================================================
if want figures; then
  log "figures: TLR5 figure inputs in figures/tlr5/"
  bash "$S/tlr5_figure_inputs.sh"
  cp macse/out/tlr5_AA.fasta figures/tlr5/
  cp "$S/tlr5_main_figure.R" "$S/tlr5_figure.R" figures/tlr5/
  ls figures/tlr5
fi


# =============================================================================
# Section 11: tlr1
# =============================================================================
if want tlr1; then
  log "tlr1: TLR1-family cluster synteny, figures/tlr1/"
  bash "$S/tlr1_cluster_inputs.sh"
  mkdir -p figures/tlr1
  cp tlr1_synteny/tlr1_synteny.tsv "$S/tlr1_synteny_figure.R" "$S/tlr1_synteny_scaled.R" figures/tlr1/
  ls figures/tlr1
fi


# =============================================================================
# Section 12: draw
# =============================================================================
if want draw; then
  log "draw: R figures"
  if command -v Rscript >/dev/null && Rscript -e \
       'quit(status = as.integer(!all(sapply(c("ggplot2","patchwork","gggenes"), requireNamespace, quietly = TRUE))))' \
       >/dev/null 2>&1; then
    ( cd figures/tlr5 && Rscript tlr5_main_figure.R && Rscript tlr5_figure.R )
    ( cd figures/tlr1 && Rscript tlr1_synteny_figure.R && Rscript tlr1_synteny_scaled.R )
  else
    echo "R with ggplot2, patchwork and gggenes is not available here. On the laptop:"
    echo "  scp -r <user>@<cluster>:$BASE/figures/tlr5 <user>@<cluster>:$BASE/figures/tlr1 <local folder>/"
    echo "  cd <local folder>/tlr5 && Rscript tlr5_main_figure.R && Rscript tlr5_figure.R"
    echo "  cd <local folder>/tlr1 && Rscript tlr1_synteny_figure.R && Rscript tlr1_synteny_scaled.R"
  fi
fi


# =============================================================================
# Section 13: record
# =============================================================================
if want record; then
  log "record: versions, environments, checksums"
  R=$REC/record_$(date +%Y%m%d_%H%M).txt
  {
    echo "TLR pipeline record: $(date)  host: $(hostname)"
    echo
    echo "== tool versions"
    conda activate miniprot
    echo "miniprot  $(miniprot --version 2>&1)"
    echo "mafft     $(mafft --version 2>&1 | head -1)"
    echo "samtools  $(samtools --version 2>&1 | head -1)"
    conda deactivate
    conda activate exonerate_env
    echo "exonerate $(exonerate --version 2>&1 | head -1)"
    conda deactivate
    echo "macse     $(conda list -n macse 2>/dev/null | awk '$1=="macse"{print $2}')"
    echo "agat      module agat/1.2.0"
    echo "python    $(python3 --version 2>&1)"
    command -v Rscript >/dev/null && echo "R         $(Rscript --version 2>&1 | head -1)"
    echo
    echo "== input checksums"
    md5sum gaga_prot.faa gaga_genomic.gff refs/chicken/cds_from_genomic.fna \
           "$TAGU/protein.faa" "$TAGU/cds_from_genomic.fna" 2>/dev/null
    echo
    echo "== output checksums"
    md5sum tlr5/results/* macse/out/tlr5_NT.fasta macse/out/tlr5_AA.fasta \
           macse/out/tlr5_disruptions.txt tlr_table/*.tsv figures/tlr5/tlr5_synteny.tsv \
           tlr1_synteny/tlr1_synteny.tsv annotation/tlr_synteny/all_tlr_flanks.tsv 2>/dev/null
  } > "$R"
  mkdir -p "$REC/envs"
  for e in miniprot exonerate_env macse; do
    conda env export -n "$e" > "$REC/envs/$e.yml" 2>/dev/null || true
  done
  echo "wrote $R and $REC/envs/"
fi


# =============================================================================
# Section 14: summary
# =============================================================================
log "summary"
[ -s tlr_table/tlr_identity_similarity.tsv ] && { echo "TLR complement (identity / similarity %):"; column -t tlr_table/tlr_identity_similarity.tsv | cut -c1-140; }
[ -s tlr5/results/tlr5_disruptions.tsv ] && { echo; echo "TLR5 disruptions (miniprot, chicken query):"; column -t tlr5/results/tlr5_disruptions.tsv; }
[ -s tlr5/results/tlr5_stops.txt ] && { echo; grep -A10 "premature stops shared" tlr5/results/tlr5_stops.txt; }
[ -s macse/out/tlr5_disruptions.txt ] && { echo; grep -A15 "^Shared by all" macse/out/tlr5_disruptions.txt; }
echo; echo "done: $(date)"