# Step 12a: Repeat curation (MCHelper v1.7.1)

git clone https://github.com/gonzalezlab/MCHelper.git
cd MCHelper
docker build -t mchelper .

docker run -it \
  -v "$(pwd)":/work \
  mchelper /bin/bash

# Download and prepare HMM profiles
wget https://urgi.versailles.inrae.fr/download/repet/profiles/ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz
tar xvf ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz
mv ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm Pfam35.0.hmm

# Combine and press BUSCO HMMs
cat hmms/*.hmm > passeriformes_odb12.hmm
hmmpress passeriformes_odb12.hmm

# Run MCHelper curation
python3 /opt/mchelper/MCHelper.py \
  -l legoRepeats-families.fa \
  -g ragtag.scaffold.Hap1.fasta \
  -o mchelp_out \
  -t 16 \
  -e 1000 \
  --input_type fasta \
  -a F \
  -r A \
  -b passeriformes_odb12.hmm

# terminal output: 12_repeatCurationMCHelper.txt
# Step 13a: Soft-masking repeats (repeatmasker v4.1.8)

bash reformatHeaders.sh

# Mask with MCHelper curated library
RepeatMasker -pa 4 -a -gff -dir repeatMaskerOutputlego -lib curated_sequences_NR_formatted.fa ragtag.scaffold.Hap1.fasta

# Mask with RepeatModeler2 raw library
RepeatMasker -pa 6 -a -gff -dir repeatMaskerOutputlego -lib legoRepeats-families.fa ragtag.scaffold.Hap1.fasta

# Compare MCHelper vs RM2 overlap
cut -f1,4,5 ragtag.scaffold.Hap1.fasta.MCH.out.gff > MCH.bed
cut -f1,4,5 ragtag.scaffold.Hap1.fasta.RM2.out.gff > RM2.bed

sort -k1,1 -k2,2n MCH.bed > MCH.sorted.bed
sort -k1,1 -k2,2n RM2.bed > RM2.sorted.bed

bedtools intersect -a MCH.sorted.bed -b RM2.sorted.bed > shared.bed
bedtools jaccard -a MCH.sorted.bed -b RM2.sorted.bed

##################################################################
# Step 13b: Dating complete repeat families

export PERL5LIB=/home/triciavanlaar/miniforge3/envs/genomics_lr/share/RepeatMasker

# All repeats (including incomplete)
perl /home/triciavanlaar/miniforge3/envs/genomics_lr/bin/calcDivergenceFromAlign.pl \
  -s legoScaffoldsHap1k2kIncomplete.divsum \
  repeatMaskerOutputlego/ragtag.scaffold.Hap1.fasta.align

# Complete repeats only (exclude _inc and _unconfirmed)
grep -v "_inc\|_unconfirmed" ragtag.scaffold.Hap1.fasta.align > legoScaffoldsHap1.complete.align

perl /home/triciavanlaar/miniforge3/envs/genomics_lr/bin/calcDivergenceFromAlign.pl \
  -s legoScaffoldsHap1k2k.complete.divsum \
  repeatMaskerOutputlego/legoScaffoldsHap1.complete.align

# Generate repeat landscape HTML
createRepeatLandscape.pl -g 1253257106 -div legoScaffoldsHap1k2kIncomplete.divsum > repeatLandscape.html

##################################################################
# R: Plot repeat landscapes

# setwd("/home/triciavanlaar/kimuraTest")

# --- All repeats landscape ---
# install.packages(c("tidyverse", "reshape2", "viridis"))

# library(tidyverse); library(reshape2); library(viridis)
# genome_size <- 1253257106
# landscape <- read.table("kimuraLandscape.tbl", header=TRUE, sep="", check.names=FALSE)
# landscape_long <- melt(landscape, id.vars="Div", variable.name="Class", value.name="Bases")
# landscape_long$BroadClass <- case_when(
#   grepl("LINE/CR1|LINE/R2|LINE$", landscape_long$Class, ignore.case=TRUE) ~ "LINE",
#   grepl("LTR/GYPSY|LTR/ERV|LTR/LARD|LTR$", landscape_long$Class, ignore.case=TRUE) ~ "LTR",
#   grepl("SINE|SINE$", landscape_long$Class, ignore.case=TRUE) ~ "SINE",
#   grepl("TIR/MERLIN|TIR/HAT|HELITRON|MAVERICK|TIR$|CLASSII", landscape_long$Class, ignore.case=TRUE) ~ "DNA",
#   grepl("Simple", landscape_long$Class, ignore.case=TRUE) ~ "Simple",
#   grepl("Unknown|#Unknown/Unknown", landscape_long$Class, ignore.case=TRUE) ~ "Unknown",
#   TRUE ~ "Other"
# )
# landscape_summary <- landscape_long %>%
#   group_by(Div, BroadClass) %>%
#   summarise(Bases=sum(Bases), .groups="drop") %>%
#   mutate(PercentGenome=Bases/genome_size*100)
# ggplot(landscape_summary, aes(x=Div, y=PercentGenome, fill=BroadClass)) +
#   geom_bar(stat="identity", position="stack", color="black", size=0.2) +
#   scale_fill_viridis(discrete=TRUE, option="D") +
#   labs(title="Repeat Landscape", x="Kimura Substitution Level (%)",
#        y="Percent of Genome Masked", fill="Repeat Class") +
#   coord_cartesian(xlim=c(0,50)) + theme_classic()
# ggsave("repeatLandscape.png", width=12, height=4, dpi=300)

# --- MCHelper complete repeats landscape ---
# filtered <- read.table("kimuraFiltered.tbl", header=TRUE, sep="", check.names=FALSE)
# filtered_long <- melt(filtered, id.vars="Div", variable.name="Class", value.name="Bases") %>%
#   mutate(SuperClass=gsub("/.*","",Class),
#     BroadClass=case_when(
#       SuperClass %in% c("LINE") ~ "LINE", SuperClass %in% c("LTR","LARD") ~ "LTR",
#       SuperClass %in% c("SINE") ~ "SINE",
#       SuperClass %in% c("TIR","HELITRON","MAVERICK","HAT","MERLIN","CLASSII") ~ "DNA",
#       grepl("Simple",Class,ignore.case=TRUE) ~ "Simple",
#       SuperClass %in% c("Unknown") ~ "Unknown", TRUE ~ "Other"))
# filtered_summary <- filtered_long %>%
#   group_by(Div, BroadClass) %>%
#   summarise(Bases=sum(Bases),.groups="drop") %>%
#   mutate(PercentGenome=Bases/genome_size*100)
# ggplot(filtered_summary, aes(x=Div,y=PercentGenome,fill=BroadClass)) + ...
# ggsave("repeatLandscapeFilt.png", width=12, height=4, dpi=300)

# --- RM2 vs MCHelper composition comparison ---
# repeat_data <- data.frame(
#   Method=rep(c("RepeatModeler2","MCHelper"),each=6),
#   Class=rep(c("LTR","LINE","SINE","DNA","Other","Unknown"),times=2),
#   Percent=c(4.53,3.71,0.01,0.03,2.94,12.00, 14.74,3.21,0.01,0.00,2.74,2.21)) %>%
#   mutate(Method=fct_rev(fct_relevel(Method,"RepeatModeler2","MCHelper")),
#          Class=factor(Class,levels=c("LTR","LINE","SINE","DNA","Other","Unknown")))
# ggplot(repeat_data,aes(x=Method,y=Percent,fill=Class)) +
#   geom_bar(stat="identity") + coord_flip() + scale_fill_brewer(palette="Set2")
# ggsave("repeatCompositionHorizontal.png", width=12, height=4, dpi=300)
