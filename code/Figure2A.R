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
ps_rar <- rarefy_even_depth(filtered_physeq, sample.size = 500, rngseed = 1337)
ps_rar 

sort(sample_sums(ps_rar))
### Subset the data for samples_types ###
food_through_time <- subset_samples(filtered_physeq, sample_type %in% c("Food") 
                                    & Nest_id %in% c("5", "13", "14", "27", "15")
                                    & Cell_ID %in% c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A"))
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "5")] <- "1"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "13")] <- "2"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "14")] <- "3"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "27")] <- "4"

# Obtain top 20 genera ####
ps_Genus2 <- tax_glom(food_through_time, taxrank = "Genus", NArm = FALSE)

# Ensure the taxonomy table has a Family column
if ("Family" %in% colnames(tax_table(ps_Genus2))) {
  
  # Replace NA values in the Genus column with "Family_unknown"
  tax_table(ps_Genus2)[is.na(tax_table(ps_Genus2)[, "Genus"]), "Genus"] <- 
    paste(as.character(tax_table(ps_Genus2)[is.na(tax_table(ps_Genus2)[, "Genus"]), "Family"]), "unclassified", sep = "_")
  
} else {
  # If the Family column doesn't exist, just replace NA with "Unknown"
  tax_table(ps_Genus2)[is.na(tax_table(ps_Genus2)[, "Genus"]), "Genus"] <- "Unclassified"
}

top20Genus2 = names(sort(taxa_sums(ps_Genus2), TRUE)[1:20])
taxtab20 = cbind(tax_table(ps_Genus2), Genus_20 = NA)
taxtab20[top20Genus2, "Genus_20"] <- as(tax_table(ps_Genus2)
                                       [top20Genus2, "Genus"], "character")
tax_table(ps_Genus2) <- tax_table(taxtab20)
ps_Genus2_ra <- transform_sample_counts(ps_Genus2, function(x) 100 * x/sum(x))
df_Genus2 <- psmelt(ps_Genus2_ra)
df_Genus2 <- arrange(df_Genus2, sample_type)
df_Genus2$Genus_20[is.na(df_Genus2$Genus_20)] <- c("Other")

# View the result
print(unique(df_Genus2$Genus_20))

# % of reads that make up the top 20 genera
mean(sample_sums(prune_taxa(top20Genus2, ps_Genus2_ra)))


### Relative abundance ####
# Custom order for Cell_ID from youngest to oldest
custom_Cell_ID_order <- c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A")
df_Genus2$Cell_ID <- factor(df_Genus2$Cell_ID, levels = custom_Cell_ID_order)

# Plot the relative abundance for pollen provisions
(fig2 <- df_Genus2 %>%
  mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my_palette) +
  facet_nested(~ Nest_id + Year + Cell_ID + sample_type + sampleid, scales = "free", space = "free") +
  labs(x = "sample_type", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size=14, face = "bold"),
    axis.title.y = element_text(size=14, face = "bold"),
    axis.ticks.y = element_line(linewidth = 1),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_textbox_simple(
      padding = margin(5, 0, 5, 0),
      margin = margin(5, 5, 5, 5),
      size = 10,
      halign = 0.5,
      fill = "white",
      box.color = "grey",
      linewidth = 1.5,
      linetype = "solid"),
    panel.background = element_blank()
  ))

ggsave("figures/legend.png", height=15, width=15)
