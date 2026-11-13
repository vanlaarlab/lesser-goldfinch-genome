#!/usr/bin/env bash
##################################################################################################################################
# Step 14: Repeat masking with and without sex chromosomes (minimap2 v2.24; RepeatMasker v4.1.5)
##################################################################################################################################
# Supporting files produced (LeGo_supporting_files/):
#   09_sexchr_scaffold_assignment.tsv              Z/W assignment of every scaffold
#   10_sexchr_assignment_summary.txt               counts and lengths per class (22 Z, 5 W)
#   11_sexchr_removed_scaffolds.txt                the 27 scaffolds excluded as Z or W
#   14_repeats_library_provenance.tsv              source and checksum of the library used for masking
#   16_repeats_repeatmasker_hits.out.gz            all RepeatMasker hits
#   17_repeats_repeatmasker_summary.tbl            RepeatMasker summary (31.71% masked)
#   18_repeats_content_by_chromosome_category.tsv  repeat content by chromosome category (Table S2)
#   19_repeats_content_with_and_without_ZW.tsv     repeat content with and without Z/W
# Input:  genomes/lego/lego_pri.fna (+ .fai) and genomes/hofi/hofi_pri.fna (step 7);
#         te_lego/library/lego_curated_sequences_NR.fa (curated library from step 13, as deposited on Dryad)
# Output: sexchr/ (Z/W assignment) and te_lego/repeatmasker/ (RepeatMasker output, summaries with and without Z/W)
# Submit: BASE=$PWD sbatch --cpus-per-task=32 --mem=64G --time=2-00:00:00 scripts/14_repeat_mask_sex_chromosomes.sh
# Conda environments: liftoff (provides minimap2 2.24) and repeatmasker (RepeatMasker 4.1.5, bioconda).

set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
cd "$BASE"
module load miniconda3/3.13
eval "$(conda shell.bash hook)"
mkdir -p sexchr te_lego/repeatmasker

GENOME=genomes/lego/lego_pri.fna
REF=genomes/hofi/hofi_pri.fna                          # house finch: CM050264.1 = Z, CM050263.1 = W
LIB=te_lego/library/lego_curated_sequences_NR.fa       # curated lego library (Dryad; https://datadryad.org/downloads/file_stream/4681394)

# 1. Align the lego scaffolds to the house finch assembly
conda activate liftoff
minimap2 --version
grep '^>' "$REF" > sexchr/hofi_pri.headers.txt
minimap2 -x asm20 --secondary=no -t "$SLURM_CPUS_PER_TASK" "$REF" "$GENOME" > sexchr/lego_pri.vs_hofi.paf
conda deactivate

# 2. Assign scaffolds to Z, W, autosomal, unplaced, mixed or unaligned
#    (MAPQ >= 5; Z/W when >= 50% of aligned bases hit Z+W; unaligned if < 5 kb or < 5% aligned)
cd sexchr
python3 - lego_pri.vs_hofi.paf hofi_pri.headers.txt "$BASE/$GENOME" lego_pri 5 0.5 5000 0.05 << 'PY'
import re, sys
from collections import defaultdict

paf, refh, qfa, name = sys.argv[1:5]
MIN_MAPQ, MIN_SHARE = int(sys.argv[5]), float(sys.argv[6])
MIN_ALN_BP, MIN_ALN_FRAC = int(sys.argv[7]), float(sys.argv[8])

# reference sequence categories from the GenBank headers
cat = {}
for line in open(refh):
    tid = line[1:].split()[0]
    m = re.search(r"chromosome (\S+?)(?:,|\s|$)", line)
    if m:
        c = m.group(1)
        cat[tid] = c if c in ("Z", "W") else "autosome"
    else:
        cat[tid] = "unplaced"

def union_bp(iv):
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if cs is None or s > ce:
            if cs is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if cs is not None else 0)

qlen = {}
q_iv = defaultdict(lambda: defaultdict(list))     # query -> category -> intervals
q_all = defaultdict(list)                          # query -> all intervals
q_tbp = defaultdict(lambda: defaultdict(int))      # query -> target -> bp
t_iv = defaultdict(list)                           # Z/W target -> intervals
n_used = n_skip = 0
for line in open(paf):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        continue
    if int(f[11]) < MIN_MAPQ or "tp:A:S" in f[12:]:
        n_skip += 1
        continue
    n_used += 1
    q, ql, qs, qe, t = f[0], int(f[1]), int(f[2]), int(f[3]), f[5]
    c = cat.get(t, "unplaced")
    qlen[q] = ql
    q_iv[q][c].append((qs, qe))
    q_all[q].append((qs, qe))
    q_tbp[q][t] += qe - qs
    if c in ("Z", "W"):
        t_iv[t].append((int(f[7]), int(f[8])))

def classify(q, length):
    aligned = union_bp(q_all[q]) if q in q_all else 0
    bp = {c: union_bp(q_iv[q][c]) for c in ("Z", "W", "autosome", "unplaced")} if q in q_iv \
         else {c: 0 for c in ("Z", "W", "autosome", "unplaced")}
    top = max(q_tbp[q].items(), key=lambda x: x[1])[0] if q in q_tbp else "NA"
    if aligned < MIN_ALN_BP or aligned < MIN_ALN_FRAC * length:
        klass = "unaligned"
    elif (bp["Z"] + bp["W"]) >= MIN_SHARE * aligned:
        klass = "Z" if bp["Z"] >= bp["W"] else "W"
    elif bp["autosome"] >= MIN_SHARE * aligned:
        klass = "autosomal"
    elif bp["unplaced"] >= MIN_SHARE * aligned:
        klass = "unplaced"
    else:
        klass = "mixed"
    return aligned, bp, top, klass

keep_fh = open(f"{name}.noZW.fna", "w")
sex_fh = open(f"{name}.ZW.fna", "w")
tab = open(f"{name}.sexchr_assignment.tsv", "w")
tab.write("scaffold\tlength\taligned_bp\taligned_frac\tZ_bp\tW_bp\tautosome_bp\t"
          "unplaced_bp\ttop_hofi_target\ttop_target_category\tclass\tremoved\n")
totals = defaultdict(lambda: [0, 0])

def finish(q, length):
    aligned, bp, top, klass = classify(q, length)
    removed = klass in ("Z", "W")
    tab.write("\t".join(map(str, [
        q, length, aligned, f"{aligned / length:.3f}" if length else "0",
        bp["Z"], bp["W"], bp["autosome"], bp["unplaced"], top,
        cat.get(top, "NA"), klass, "yes" if removed else "no"])) + "\n")
    totals[klass][0] += 1
    totals[klass][1] += length

cur, cur_len, out = None, 0, None
with open(qfa) as fh:
    for line in fh:
        if line.startswith(">"):
            if cur is not None:
                finish(cur, cur_len)
            cur, cur_len = line[1:].split()[0], 0
            # the class only needs the PAF and the length; unaligned scaffolds
            # are never removed, so the output file can be chosen now
            length_hint = qlen.get(cur, 0)
            _, _, _, k = classify(cur, length_hint) if length_hint else (0, 0, 0, "unaligned")
            out = sex_fh if k in ("Z", "W") else keep_fh
        else:
            cur_len += len(line.strip())
        out.write(line)
if cur is not None:
    finish(cur, cur_len)
keep_fh.close(); sex_fh.close(); tab.close()

ref_len = {}
for line in open(paf):
    f = line.split("\t")
    if len(f) > 6:
        ref_len[f[5]] = int(f[6])
lines = [f"{name} vs house finch: {n_used} alignments used, {n_skip} skipped "
         f"(MAPQ < {MIN_MAPQ} or secondary)",
         f"rules: removed if Z+W >= {MIN_SHARE:.0%} of aligned bases; unaligned if "
         f"< {MIN_ALN_BP} bp or < {MIN_ALN_FRAC:.0%} of length aligned", "",
         f"{'class':<10} {'scaffolds':>9} {'Mb':>9}  action"]
for k in ("autosomal", "Z", "W", "unplaced", "mixed", "unaligned"):
    n, b = totals.get(k, [0, 0])
    lines.append(f"{k:<10} {n:>9} {b / 1e6:>9.2f}  {'REMOVED' if k in ('Z', 'W') else 'kept'}")
lines.append("")
for t in sorted(t_iv):
    cov = union_bp(t_iv[t])
    L = ref_len.get(t, 0)
    lines.append(f"house finch {cat[t]} ({t}, {L / 1e6:.1f} Mb): "
                 f"{cov / 1e6:.1f} Mb covered by lego alignments "
                 f"({cov / L:.0%})" if L else f"house finch {t}: no length")
summary = "\n".join(lines)
open(f"{name}.sexchr_summary.txt", "w").write(summary + "\n")
print(summary)
PY
awk -F'\t' 'NR>1 && $12=="yes"{print $1}' lego_pri.sexchr_assignment.tsv > lego_pri.removed_scaffolds.txt
cd "$BASE"

# 3. RepeatMasker on the whole assembly with the curated library
conda activate repeatmasker
unset PERL5LIB PERLLIB
RepeatMasker -v
RepeatMasker -pa $((SLURM_CPUS_PER_TASK / 4)) -lib "$LIB" -xsmall -gff -a -no_is -dir te_lego/repeatmasker "$GENOME"

# 4. Exclude the Z/W scaffolds and summarize repeat content with and without them
cp sexchr/lego_pri.sexchr_assignment.tsv te_lego/repeatmasker/
python3 - "$GENOME.fai" sexchr/lego_pri.sexchr_assignment.tsv \
          te_lego/repeatmasker/lego_pri.fna.out te_lego/repeatmasker/lego_pri.fna.out.gff \
          te_lego/repeatmasker/lego_pri.fna.masked te_lego/repeatmasker << 'PY'
import sys
from collections import defaultdict
fai, assign, rm_out, rm_gff, rm_masked, outdir = sys.argv[1:7]
MAJOR = ["LINE", "SINE", "LTR", "DNA", "RC", "Retroposon", "ClassI_unknown", "Unknown",
         "Satellite", "Simple_repeat", "Low_complexity"]
CATS = ["autosomal", "Z", "W", "other", "noZW", "all"]

length = {}
for line in open(fai):
    n, l = line.split("\t")[:2]
    length[n] = int(l)
cat = {}
for i, line in enumerate(open(assign)):
    if i == 0: continue
    f = line.rstrip("\n").split("\t")
    cat[f[0]] = {"autosomal": "autosomal", "Z": "Z", "W": "W"}.get(f[10], "other")
zw = {n for n, c in cat.items() if c in ("Z", "W")}

# library labels are Wicker style (CLASSI/LTR, CLASSII/TIR, ...); RepeatMasker's
# own simple repeats keep RepeatMasker labels. Grouped the way RepeatMasker does.
WICKER = {"LTR": "LTR", "DIRS": "LTR", "LINE": "LINE", "PLE": "LINE", "PENELOPE": "LINE",
          "SINE": "SINE", "TIR": "DNA", "MITE": "DNA", "CRYPTON": "DNA",
          "MAVERICK": "DNA", "HELITRON": "RC"}
other_labels = defaultdict(int)
def to_major(label):
    parts = label.split("/")
    if parts[0].upper().startswith("CLASS"):
        sub = parts[1].upper() if len(parts) > 1 else ""
        m = WICKER.get(sub)
        if m is None and sub in ("", "UNKNOWN", "UNCLASSIFIED", "NOCAT"):
            m = "ClassI_unknown" if parts[0].upper() == "CLASSI" else "DNA"
    else:
        m = parts[0] if parts[0] in MAJOR else None
    if m is None:
        other_labels[label] += 1
        return "Other"
    return m

def union_bp(iv):
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if cs is None or s > ce:
            if cs is not None: tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if cs is not None else 0)

iv_all = defaultdict(list)
iv_cls = defaultdict(lambda: defaultdict(list))
with open(rm_out) as fin, open(f"{outdir}/lego_pri.noZW.fna.out", "w") as fout:
    for i, line in enumerate(fin):
        f = line.split()
        if i < 3:
            fout.write(line); continue
        if len(f) < 11 or not f[0].isdigit(): continue
        q, s, e = f[4], int(f[5]) - 1, int(f[6])
        if q not in zw: fout.write(line)
        major = to_major(f[10])
        iv_all[q].append((s, e)); iv_cls[q][major].append((s, e))
with open(rm_gff) as fin, open(f"{outdir}/lego_pri.noZW.fna.gff", "w") as fout:
    for line in fin:
        if line.startswith("#") or line.split("\t", 1)[0] not in zw:
            fout.write(line)
with open(rm_masked) as fin, open(f"{outdir}/lego_pri.noZW.masked.fna", "w") as fout:
    keep = True
    for line in fin:
        if line.startswith(">"):
            keep = line[1:].split()[0] not in zw
        if keep: fout.write(line)

acc = {c: defaultdict(float) for c in CATS}
for n, L in length.items():
    c = cat.get(n, "other")
    groups = [c, "all"] + (["noZW"] if c in ("autosomal", "other") else [])
    masked = union_bp(iv_all[n]) if n in iv_all else 0
    per = {k: union_bp(v) for k, v in iv_cls[n].items()} if n in iv_cls else {}
    for g in groups:
        acc[g]["n"] += 1; acc[g]["len"] += L; acc[g]["masked"] += masked
        for k, v in per.items(): acc[g][k] += v

head = ["assembly", "category", "scaffolds", "length_Mb", "masked_Mb", "masked_pct"] + \
       [f"{k}_pct" for k in MAJOR + ["Other"]]
rows = []
for c in CATS:
    a = acc[c]
    if not a["n"]: continue
    rows.append(["lego_pri", c, int(a["n"]), f"{a['len'] / 1e6:.2f}", f"{a['masked'] / 1e6:.2f}",
                 f"{100 * a['masked'] / a['len']:.2f}"] +
                [f"{100 * a[k] / a['len']:.2f}" for k in MAJOR + ["Other"]])
with open(f"{outdir}/te_summary_by_category.tsv", "w") as out:
    out.write("\t".join(head) + "\n")
    for r in rows: out.write("\t".join(map(str, r)) + "\n")

sens = []
for k in ["total"] + MAJOR + ["Other"]:
    key = "masked" if k == "total" else k
    all_pct = 100 * acc["all"][key] / acc["all"]["len"]
    nzw_pct = 100 * acc["noZW"][key] / acc["noZW"]["len"] if acc["noZW"]["len"] else 0.0
    share = 100 * (acc["Z"][key] + acc["W"][key]) / acc["all"][key] if acc["all"][key] else 0.0
    sens.append(["lego_pri", k, f"{all_pct:.2f}", f"{nzw_pct:.2f}", f"{nzw_pct - all_pct:+.2f}", f"{share:.1f}"])
with open(f"{outdir}/te_sensitivity_noZW.tsv", "w") as out:
    out.write("assembly\tclass\tall_pct\tnoZW_pct\tdifference_pp\tshare_on_ZW_pct\n")
    for r in sens: out.write("\t".join(r) + "\n")

if other_labels:
    print("labels counted as 'Other': " + ", ".join(
        f"{k} ({v})" for k, v in sorted(other_labels.items(), key=lambda x: -x[1])[:20]))
print(f"excluded {len(zw)} Z/W scaffolds ({sum(length.get(n, 0) for n in zw) / 1e6:.1f} Mb)")
print(f"\n{'category':<10} {'scaffolds':>9} {'Mb':>9} {'masked%':>8} {'LINE%':>6} {'LTR%':>6} {'DNA%':>6}")
for r in rows:
    d = dict(zip(head, r))
    print(f"{d['category']:<10} {d['scaffolds']:>9} {d['length_Mb']:>9} {d['masked_pct']:>8} "
          f"{d['LINE_pct']:>6} {d['LTR_pct']:>6} {d['DNA_pct']:>6}")
print("\nwith and without Z/W (percent of assembly; pp = percentage points)")
print(f"{'class':<15} {'all%':>7} {'noZW%':>7} {'diff pp':>8} {'on Z+W %':>9}")
for r in sens:
    if r[1] in ("total", "LINE", "SINE", "LTR", "DNA", "RC", "ClassI_unknown", "Satellite", "Simple_repeat"):
        print(f"{r[1]:<15} {r[2]:>7} {r[3]:>7} {r[4]:>8} {r[5]:>9}")
PY