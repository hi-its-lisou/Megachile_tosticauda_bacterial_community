### Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(extrafont)
library(forcats)
library(ggforce)
library(vegan)
library(xfun)

### Preparing the data ###
# Load phyloseq object
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_5.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps #1089 taxa and 208 samples

# Distinguish the two most prevalent ASVs
tax_table_df <- as.data.frame(tax_table(ps))
# Rename the genus using *
tax_table_df["c624f2e4228eea7296b2a77e2d4b7e50", "Genus"] <- "Acinetobacter*"
tax_table_df["80626c0d45293428d118ce1f05a1ab18", "Genus"] <- "Tyzzerella*"
# Update the taxonomy
tax_table(ps) <- as.matrix(tax_table_df)

# Load custom colour palette for each taxa
colours_df <- read_csv("input_files/colour_list.csv")
my_palette <- colours_df$colours
names(my_palette) <- colours_df$genera
my_palette

# Remove contaminants, chloroplasts, and mitochondria ####
contam <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                     delim = "\n", 
                     col_names = "asv")

chloro_mito_decontam_asvs <- contam$asv
all_asvs <- taxa_names(ps)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
filtered_physeq <- prune_taxa(asvs_to_keep, ps)
filtered_physeq

sort(sample_sums(filtered_physeq))


# Rarefy the data
ps_rar <- rarefy_even_depth(filtered_physeq, sample.size = 500, rngseed = 1337) #1500 as determined by QIIME output and plateau of curves
ps_rar 

sort(sample_sums(ps_rar))

# Subset the data for samples_types ####
food_with_brood_through_time <- subset_samples(filtered_physeq, sample_type %in% c("Food") 
                                    & Nest_id %in% c("5", "13", "14", "27")
                                    & Cell_ID %in% c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A"))
food_with_brood_through_time@sam_data$Nest_id[which(food_with_brood_through_time@sam_data$Nest_id == "5")] <- "1"
food_with_brood_through_time@sam_data$Nest_id[which(food_with_brood_through_time@sam_data$Nest_id == "13")] <- "2"
food_with_brood_through_time@sam_data$Nest_id[which(food_with_brood_through_time@sam_data$Nest_id == "14")] <- "3"
food_with_brood_through_time@sam_data$Nest_id[which(food_with_brood_through_time@sam_data$Nest_id == "27")] <- "4"



# Alpha diversity of pollen provisions as they are consumed ####
alpha_diversity <- alpha(food_with_brood_through_time, index = "Shannon")
metadata <- meta(food_with_brood_through_time)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadata <- merge(alpha_diversity, metadata, by = "name")

# Set the factor levels for Cell_ID in the merged dataframe
custom_Cell_ID_order <- c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A")
alpha_diversity_metadata$Cell_ID <- factor(alpha_diversity_metadata$Cell_ID, levels = custom_Cell_ID_order)

# Create the plot with regression line and R-squared annotation
boxplot <- alpha_diversity_metadata %>%
  ggplot(aes(x = Cell_ID, y = diversity_shannon, colour = Cell_ID)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(y = "Shannon diversity", x = "Cell position (youngest to oldest)") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = 'bold'),
        legend.position = "none")

# Convert levels to numeric values for regression
numeric_cell_ids <- as.numeric(factor(alpha_diversity_metadata$Cell_ID, levels = custom_Cell_ID_order))

# Create the plot
ggplot(alpha_diversity_metadata, aes(x = numeric_cell_ids, y = diversity_shannon)) +
  geom_point(size = 5, color = "black") +  
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  labs(y = "Shannon diversity", x = "Cell position (youngest to oldest)") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = 'bold')) +
  facet_wrap(~ Nest_id) +
  scale_x_continuous(breaks = 1:length(custom_Cell_ID_order), 
                     labels = custom_Cell_ID_order)


ggsave("figures/2B.png", height=8, width=12)


model_interaction <- lm(diversity_shannon ~ as.numeric(Cell_ID) * Nest_id, data = alpha_diversity_metadata)
summary(model_interaction)
