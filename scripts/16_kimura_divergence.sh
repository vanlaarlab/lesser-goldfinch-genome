#!/bin/bash
# Step 16: Kimura divergence for the repeat landscapes (RepeatMasker calcDivergenceFromAlign.pl),
# for the full lesser goldfinch primary assembly and without the Z/W scaffolds (step 14).
#
# Usage, from the project directory:   BASE=$PWD bash scripts/16_kimura_divergence.sh
# Outputs in te_lego/repeatmasker/:
#   lego_pri.divsum, lego_pri.noZW.divsum             full calcDivergenceFromAlign summaries
#   lego_pri.kimura.tbl, lego_pri.noZW.kimura.tbl     landscape tables (Div x repeat class)
# Supporting files produced (LeGo_supporting_files/):
#   20_repeats_divergence_summary_full_assembly.divsum  divergence summary, full assembly
#   21_repeats_divergence_summary_without_ZW.divsum     divergence summary, without Z/W
#   22_repeats_landscape_full_assembly.tsv              landscape table, full assembly
#   23_repeats_landscape_without_ZW.tsv                 landscape table, without Z/W (Figure 1)

set -euo pipefail
BASE=${BASE:?set BASE to the project directory}
cd "$BASE/te_lego/repeatmasker"
module load miniconda3/3.13
eval "$(conda shell.bash hook)"
conda activate repeatmasker
unset PERL5LIB PERLLIB

# Use the copy inside RepeatMasker's own folder, next to its Perl modules
# (the copy in the env's bin/ cannot find them)
RMDIR=$(dirname "$(find "$CONDA_PREFIX" -name SearchResult.pm -print -quit)")
CALC=$RMDIR/util/calcDivergenceFromAlign.pl
[ -s "$RMDIR/SearchResult.pm" ] && [ -s "$CALC" ] \
  || { echo "RepeatMasker util scripts not found under $CONDA_PREFIX" >&2; exit 1; }
export PERL5LIB=$RMDIR

# 1. Full assembly
perl "$CALC" -s lego_pri.divsum lego_pri.fna.align

# 2. Without Z/W: drop every alignment whose query is a scaffold removed in Step 18.
#    Each alignment in a RepeatMasker .align file starts with a header line
#    (score, %div, %del, %ins, query, begin, end, ...); field 5 is the query.
awk 'NR == FNR { zw[$1] = 1; next }
     /^ *[0-9]+ +[0-9.]+ +[0-9.]+ +[0-9.]+ +[^ ]+ +[0-9]+ +[0-9]+ / { keep = !($5 in zw) }
     keep' "$BASE/sexchr/lego_pri.removed_scaffolds.txt" lego_pri.fna.align > lego_pri.noZW.align
perl "$CALC" -s lego_pri.noZW.divsum lego_pri.noZW.align

# 3. Landscape tables: the Kimura section at the end of each .divsum
for f in lego_pri lego_pri.noZW; do
  sed -n '/^Coverage for each repeat class and divergence (Kimura)/,$p' "$f.divsum" | tail -n +2 > "$f.kimura.tbl"
  [ -s "$f.kimura.tbl" ] || echo "WARNING: no Kimura table found in $f.divsum" >&2
  echo "$f.kimura.tbl: $(($(wc -l < "$f.kimura.tbl") - 1)) divergence bins, $(head -1 "$f.kimura.tbl" | wc -w) columns"
done
