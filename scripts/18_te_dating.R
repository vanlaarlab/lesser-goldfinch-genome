# Step 18: Age of the most recent TE activity peak (Table 2) and the landscape values quoted in the text
#
# Run inside te_lego/repeatmasker/ after steps 14 and 16:
#   cd $BASE/te_lego/repeatmasker && Rscript $BASE/scripts/18_te_dating.R
# Inputs: lego_pri.noZW.kimura.tbl (step 16), te_summary_by_category.tsv (step 14)
#
# calcDivergenceFromAlign.pl truncates divergence to whole percents, so bin d holds copies with d to d+1%
# divergence. Ages use K = 2*mu*t (t = K / (2*mu)) at both edges of the peak bin, with per-generation rates from
# three songbird pedigree studies converted with the generation length of S. psaltria.
# Supporting files used (LeGo_supporting_files/): 18_repeats_content_by_chromosome_category.tsv, 23_repeats_landscape_without_ZW.tsv
#   output: Table 2 and the landscape values in the text

rates_per_generation <- c(
  "Collared flycatcher (Smeds et al., 2016)"           = 4.6e-9,
  "Zebra finch (Prentout et al., 2025)"                = 5.0e-9,
  "Great reed warbler, autosomes (Zhang et al., 2023)" = 7.16e-9
)
generation_length <- 2.8   # years; BirdLife International, S. psaltria

land <- read.table("lego_pri.noZW.kimura.tbl", header = TRUE, sep = "", check.names = FALSE)
summ <- read.delim("te_summary_by_category.tsv", check.names = FALSE)
assembly_bp <- summ$length_Mb[summ$category == "noZW"] * 1e6

pct <- rowSums(land[, -1]) / assembly_bp * 100
names(pct) <- land$Div
young <- land$Div[land$Div <= 5]
peak  <- young[which.max(pct[as.character(young)])]
older <- land$Div[land$Div >= 6 & land$Div <= 15]
peak2 <- older[which.max(pct[as.character(older)])]

mu  <- rates_per_generation / generation_length
age <- function(K) K / (2 * mu) / 1e6

cat(sprintf("Most recent peak: %d-%d%% bin (%.2f%% of the assembly without Z/W)\n", peak, peak + 1, pct[as.character(peak)]))
cat(sprintf("Copies within 2%% of their consensus (bins 0 and 1): %.1f%% of the assembly\n", sum(pct[c("0", "1")])))
cat(sprintf("Older LTR peak: %d-%d%% bin\n\n", peak2, peak2 + 1))

table2 <- data.frame(
  rate_source    = names(rates_per_generation),
  per_generation = rates_per_generation,
  per_year       = signif(mu, 3),
  age_My         = sprintf("%.1f-%.1f", age(peak / 100), age((peak + 1) / 100)),
  row.names = NULL
)
print(table2, right = FALSE)
cat(sprintf("\nMost recent peak, range across rates: %.1f-%.1f million years ago\n",
            min(age(peak / 100)), max(age((peak + 1) / 100))))
cat(sprintf("Older peak, range across rates: %.0f-%.0f million years ago\n",
            min(age(peak2 / 100)), max(age((peak2 + 1) / 100))))
write.table(table2, "table2_te_dating.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
