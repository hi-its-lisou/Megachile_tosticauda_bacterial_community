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


# Rarefy the data, 
ps_rar <- rarefy_even_depth(filtered_physeq, sample.size = 500, rngseed = 1337)
ps_rar 

sort(sample_sums(ps_rar))

### Subset the data for samples_types ###
#Honey bees
HB_subset <- subset_samples(ps_rar, sample_type %in% c("Honey_bee"))

#Prepupae
prepupal_subset <- subset_samples(ps_rar, sample_type %in% c("Prepupae") & is.na(AB_treatment) & !Nest_id %in% c("5", "13", "14", "27"))

#Larvae
larva_subset <- subset_samples(ps_rar, sample_type %in% c("Larvae") & Year %in% c("2022"))
mature_larvae <- subset_samples(ps_rar, sample_type %in% c("Prepupae") & Nest_id %in% c("5", "13", "14", "27"))
mature_larvae@sam_data$sample_type[which(mature_larvae@sam_data$sample_type == "Prepupae")] <- "Mature larvae"

larvae_all <- merge_phyloseq (mature_larvae, larva_subset)
larvae_all@sam_data$sample_type[which(larvae_all@sam_data$sample_type == "Mature larvae")] <- "Larvae"

#Natural adults
adults_subset <- subset_samples(ps_rar, sample_type == "Adults" & Env_exposure == "Free-flying")

#Frass
frass_subset <- subset_samples(ps_rar, Description_of_sample %in% c("CocoonFrass", "Frass", "PinkCocoonFrass"))
frass_subset@sam_data$sample_type[which(frass_subset@sam_data$sample_type == "Nest_contents")] <- "Frass contents"




### Analysis development stages ###
# Compare the beta diversity of the different development stages ####
development <- merge_phyloseq(prepupal_subset, adults_subset, larvae_all)
unweighted_unifrac <- ordinate(development, method = "PCoA", distance = "unifrac", weighted=F)

colours <- c("#f87970", "green", "#4B0092")

(Figure1B <- plot_ordination(physeq = development, 
                             ordination = unweighted_unifrac, 
                             color = "sample_type",
                             axes = c(1, 2)) +
    geom_point(size = 5, alpha = 0.6) +
    stat_ellipse(geom = "polygon", type="norm", alpha=0) +
    theme_minimal() +
    scale_color_manual(values = colours) +
    theme (axis.text.y = element_text(size=18, face = 'bold'),
           axis.title.y = element_text(size=14, face = 'bold'),
           axis.text.x = element_text(size=18, face = 'bold'),
           axis.title.x = element_text(size=14, face = 'bold'),
           legend.text = element_text(size = 14),
           legend.position = "top", 
           legend.title = element_blank()))

# Adonis statistic to determine whether communities are significantly different
#Refactor metadata
md_combined_subset <- as(sample_data(development), "data.frame")
#Run adonis with 9999 permutations
uw_adonis <- adonis2(distance(development, method="unifrac") ~ sample_type, 
                     data = md_combined_subset, 
                     permutations = 9999)
uw_adonis

#Check for homogenaeity multivariate dispersion with betadisper
homogenaeity <- md_combined_subset[["sample_type"]]
uw_disp <- betadisper(distance(development, method="unifrac"), group = homogenaeity)
anova(uw_disp)

plot(uw_disp)
boxplot(uw_disp)
#Tukey HSD
uw_disp_HSD <- TukeyHSD(uw_disp)
plot(uw_disp_HSD)


# Combine the subsetted data ####
combined_subset <- merge_phyloseq(prepupal_subset, adults_subset, larva_subset, mature_larvae, frass_subset)

# Obtain top 20 genera ####
ps_Genus1 <- tax_glom(combined_subset, taxrank = "Genus", NArm = FALSE)

if ("Family" %in% colnames(tax_table(ps_Genus1))) {
  
  # Replace NA values in the Genus column with "Family_unknown"
  tax_table(ps_Genus1)[is.na(tax_table(ps_Genus1)[, "Genus"]), "Genus"] <- 
    paste(as.character(tax_table(ps_Genus1)[is.na(tax_table(ps_Genus1)[, "Genus"]), "Family"]), "unclassified", sep = "_")
  
} else {
  # If the Family column doesn't exist, just replace NA with "Unknown"
  tax_table(ps_Genus1)[is.na(tax_table(ps_Genus1)[, "Genus"]), "Genus"] <- "Unclassified"
}

top20Genus1 = names(sort(taxa_sums(ps_Genus1), TRUE)[1:19])
taxtab20 = cbind(tax_table(ps_Genus1), Genus_20 = NA)
taxtab20[top20Genus1, "Genus_20"] <- as(tax_table(ps_Genus1)
                                       [top20Genus1, "Genus"], "character")
tax_table(ps_Genus1) <- tax_table(taxtab20)
ps_Genus1_ra <- transform_sample_counts(ps_Genus1, function(x) 100 * x/sum(x))
df_Genus1 <- psmelt(ps_Genus1_ra)
df_Genus1 <- arrange(df_Genus1, sample_type)
df_Genus1$Genus_20[is.na(df_Genus1$Genus_20)] <- c("Other")

# View the result
print(unique(df_Genus1$Genus_20))

# % of reads that make up the top 20 genera
mean(sample_sums(prune_taxa(top20Genus1, ps_Genus1_ra)))


### Relative abundance ####
# Custom order for sample types
custom_order <- c("Adults", "Larvae", "Mature_Larvae", "Prepupae", "Frass contents")
df_Genus1$sample_type <- factor(df_Genus1$sample_type, levels = custom_order)

# Custom order for cell ID from youngest to oldest
custom_Cell_ID_order <- c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A")
df_Genus1$Cell_ID <- factor(df_Genus1$Cell_ID, levels = custom_Cell_ID_order)

# Plot the relative abundance for development stages
(Figure1A <- df_Genus1 %>%
  mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my_palette) +
  facet_nested(~ sample_type + Year + Cell_ID + sampleid, 
               scales = "free", 
               space = "free") +
  labs(x = "sampleid", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size = 14, face = 'bold'),
    axis.title.y = element_text(size = 14, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    strip.text = element_textbox_simple(
      padding = margin(2, 0, 2, 0),
      size = 12,
      face = "bold",
      halign = 0.5,
      fill = "grey",
      box.color = "black",
      linewidth = 0.5,
      linetype = "solid"),
    panel.background = element_blank()
  ))


# Use cowplot to stich a/b
Figure1 <- cowplot::plot_grid(Figure1A,
                              Figure1B,
                              ncol = 2,
                              rel_heights = c(1,0.6),
                              rel_widths = c(1, 0.5))
Figure1

#Save plot
ggsave("figures/OctFigure1.png", height=10, width=18)




