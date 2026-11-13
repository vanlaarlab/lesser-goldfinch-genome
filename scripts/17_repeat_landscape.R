# Step 17: Repeat landscape and composition (Figure 1), lesser goldfinch primary assembly
#    without the Z/W scaffolds.
#
# Run inside te_lego/repeatmasker/ (or copy the inputs elsewhere):
#   cd $BASE/te_lego/repeatmasker && Rscript $BASE/scripts/17_repeat_landscape.R
# Inputs:
#   lego_pri.noZW.kimura.tbl      (kimura_divergence.sh)
#   te_summary_by_category.tsv    (Step 18)
# Outputs:
#   repeat_landscape_noZW.png / .pdf
#   repeat_composition_noZW.png / .tsv            (composition in the earlier categories)
#   repeat_landscape_and_composition_noZW.png   (both panels together)
# Supporting files used (LeGo_supporting_files/): 18_repeats_content_by_chromosome_category.tsv, 23_repeats_landscape_without_ZW.tsv
#   output: Figure 1 (repeat landscape and composition)

# install.packages(c("tidyverse", "patchwork"))
library(tidyverse)
library(patchwork)

# Same six categories, order and colors as the paper's other repeat landscape
# (viridis, option D). All six are always shown in the legend, even at 0%.
class_levels  <- c("LINE", "LTR", "Other", "Simple", "SINE", "Unknown")
class_colors <- c(LINE = "#440154", LTR = "#414487", Other = "#2A788E",
                   Simple = "#22A884", SINE = "#7AD151", Unknown = "#FDE725")

# Classification rules from the earlier repeat-landscape script; DNA elements
# (TIR, Helitron), not a category in the other figure, are counted as Other
broad_class <- function(label) {
  case_when(
    grepl("LINE/CR1|LINE/R2|LINE$", label, ignore.case = TRUE) ~ "LINE",
    grepl("LTR/GYPSY|LTR/ERV|LTR/LARD|LTR$", label, ignore.case = TRUE) ~ "LTR",
    grepl("SINE|SINE$", label, ignore.case = TRUE) ~ "SINE",
    grepl("TIR/MERLIN|TIR/HAT|HELITRON|MAVERICK|TIR$|CLASSII", label, ignore.case = TRUE) ~ "Other",
    grepl("Simple", label, ignore.case = TRUE) ~ "Simple",
    grepl("Unknown|#Unknown/Unknown", label, ignore.case = TRUE) ~ "Unknown",
    TRUE ~ "Other"
  )
}

# Length (bp) of the assembly without Z/W, from the Step 18 summary
summ <- read.delim("te_summary_by_category.tsv", check.names = FALSE)
size_nozw <- summ$length_Mb[summ$category == "noZW"] * 1e6

# --- Repeat landscape ------------------------------------------------------------
land <- read.table("lego_pri.noZW.kimura.tbl", header = TRUE, sep = "", check.names = FALSE) %>%
  pivot_longer(-Div, names_to = "Class", values_to = "Bases") %>%
  mutate(BroadClass = broad_class(Class)) %>%
  group_by(Div, BroadClass) %>%
  summarize(Bases = sum(Bases), .groups = "drop") %>%
  mutate(PercentAssembly = Bases / size_nozw * 100,
         BroadClass = factor(BroadClass, levels = class_levels))

# --- Repeat composition --------------------------------------------------------------
comp <- summ %>%
  filter(category == "noZW") %>%
  select(ends_with("_pct"), -masked_pct) %>%
  pivot_longer(everything(), names_to = "Class", values_to = "Percent") %>%
  mutate(Class = sub("_pct$", "", Class),
         Class = case_when(
           Class %in% c("LINE", "LTR", "SINE", "Unknown") ~ Class,
           Class == "Simple_repeat" ~ "Simple",
           TRUE                     ~ "Other")) %>%
  group_by(Class) %>%
  summarize(Percent = sum(Percent), .groups = "drop") %>%
  filter(Percent > 0) %>%
  mutate(Set = "Without Z/W",
         Class = factor(Class, levels = class_levels))

p_land <- ggplot(land, aes(x = Div, y = PercentAssembly, fill = BroadClass)) +
  geom_col(color = "black", linewidth = 0.2, show.legend = TRUE) +
  scale_fill_manual(values = class_colors, limits = class_levels, breaks = class_levels) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura Substitution Level (%)", y = "Percent of Genome Masked",
       fill = "Repeat Class") +
  theme_classic()

ggsave("repeat_landscape_noZW.png", p_land, width = 12, height = 4, dpi = 300, bg = "white")
ggsave("repeat_landscape_noZW.pdf", p_land, width = 12, height = 4)

p_comp <- ggplot(comp, aes(x = Set, y = Percent, fill = Class)) +
  geom_col(color = "black", linewidth = 0.2, show.legend = TRUE) +
  coord_flip() +
  scale_fill_manual(values = class_colors, limits = class_levels, breaks = class_levels) +
  labs(x = NULL, y = "Percent of assembly", fill = "Repeat Class") +
  theme_classic()

ggsave("repeat_composition_noZW.png", p_comp, width = 12, height = 2.5, dpi = 300, bg = "white")
write.table(comp %>% arrange(Class) %>% transmute(Class, Percent = round(Percent, 2)),
            "repeat_composition_noZW.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# --- Both together -----------------------------------------------------------------
p_both <- p_land / p_comp + plot_layout(heights = c(3, 1), guides = "collect")
ggsave("repeat_landscape_and_composition_noZW.png", p_both, width = 12, height = 6, dpi = 300, bg = "white")
