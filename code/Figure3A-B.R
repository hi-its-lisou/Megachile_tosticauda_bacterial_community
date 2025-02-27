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


# Subset the data for adults in the acquisition experiment ####
adults2023 <- subset_samples(ps_rar, sample_type %in% c("Adults"))
adults2023

# Obtain top 20 genera ####
ps_Genus3 <- tax_glom(adults2023, taxrank = "Genus", NArm = FALSE)

# Ensure the taxonomy table has a Family column
if ("Family" %in% colnames(tax_table(ps_Genus3))) {
  
  # Replace NA values in the Genus column with "Family_unknown"
  tax_table(ps_Genus3)[is.na(tax_table(ps_Genus3)[, "Genus"]), "Genus"] <- 
    paste(as.character(tax_table(ps_Genus3)[is.na(tax_table(ps_Genus3)[, "Genus"]), "Family"]), "unclassified", sep = "_")
  
} else {
  # If the Family column doesn't exist, just replace NA with "Unknown"
  tax_table(ps_Genus3)[is.na(tax_table(ps_Genus3)[, "Genus"]), "Genus"] <- "Unclassified"
}

top20Genus3 = names(sort(taxa_sums(ps_Genus3), TRUE)[1:19])
taxtab20 = cbind(tax_table(ps_Genus3), Genus_20 = NA)
taxtab20[top20Genus3, "Genus_20"] <- as(tax_table(ps_Genus3)
                                       [top20Genus3, "Genus"], "character")
tax_table(ps_Genus3) <- tax_table(taxtab20)
ps_Genus3_ra <- transform_sample_counts(ps_Genus3, function(x) 100 * x/sum(x))
df_Genus3 <- psmelt(ps_Genus3_ra)
df_Genus3 <- arrange(df_Genus3, sample_type)
df_Genus3$Genus_20[is.na(df_Genus3$Genus_20)] <- c("Other")

# View the result
print(unique(df_Genus3$Genus_20))

# % of reads that make up the top 20 genera ####
mean(sample_sums(prune_taxa(top20Genus3, ps_Genus3_ra)))

df_Genus3$Genus <- factor(df_Genus3$Genus, 
                         levels = c("Acinetobacter*", "Tyzzerella*", "Acinetobacter", "Tyzzerella"))

# Plot the relative abundance ####
(Figure3A <- df_Genus3 %>%
  mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my_palette) +
  facet_nested(~ sample_type + Env_exposure + Gender + sampleid.1, 
               scales = "free", space = "free") +
  labs(x = "sample_type", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size=14, face = 'bold'),
    axis.title.y = element_text(size=14, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_textbox_simple(
      padding = margin(5, 0, 5, 0),
      margin = margin(5, 5, 5, 5),
      size = 10,
      face = "bold",
      halign = 0.5,
      fill = "white",
      box.color = "grey",
      linewidth = 1.5,
      linetype = "solid"),
    panel.background = element_blank()
  ))

# Beta diversity - PCOA ordination of unweighted UniFrac distances ####
unweighted_unifrac <- ordinate(adults2023, method = "PCoA", distance = "unifrac", weighted=F)

adults_colours <- c("#f87970", "#1E88E5", "#FFC107")

(Figure3B <- plot_ordination(physeq = adults2023, 
                             ordination = unweighted_unifrac, 
                             color = "Env_exposure",
                             axes = c(1, 2)) +
    geom_point(size = 5, alpha = 0.6) +
    stat_ellipse(geom = "polygon", type="norm", alpha=0) +
        theme_minimal() +
    scale_color_manual(values = adults_colours) +
    theme (axis.text.y = element_text(size=18, face = 'bold'),
           axis.title.y = element_text(size=18, face = 'bold'),
           axis.text.x = element_text(size=18, face = 'bold'),
           axis.title.x = element_text(size=18, face = 'bold'),
           legend.text = element_text(size = 18),
           legend.position = "top", 
           legend.title = element_blank()))

# Adonis statistic to determine whether communities are significantly different
#Refactor metadata
adults_combined_subset <- as(sample_data(adults2023), "data.frame")
#Run adonis with 9999 permutations
uw_adonis <- adonis2(distance(adults2023, method="unifrac") ~ Env_exposure, 
                     data = adults_combined_subset, 
                     permutations = 9999)
uw_adonis


#Check for homogenaeity multivariate dispersion with betadisper
homogenaeity <- adults_combined_subset[["Env_exposure"]]
uw_disp <- betadisper(distance(adults2023, method="unifrac"), group = homogenaeity)
anova(uw_disp)

plot(uw_disp)
boxplot(uw_disp)
#Tukey HSD
uw_disp_HSD <- TukeyHSD(uw_disp)
plot(uw_disp_HSD)



# Combine the plots ####
AE_plot <- cowplot::plot_grid(Figure3A,
                               Figure3B,
                               ncol = 2,
                               rel_heights = c(1,0.6),
                               rel_widths = c(1, 0.5))
AE_plot
# Save plot ####
ggsave("figures/OctFigure3.png", height=13, width=20)


##For supplementary materials with read counts of Acinetobacter and Tyzzerella
# Acinetobacter
source("Plot_specific_ASV_Abundance_Function.R")
plot.specific.ASV(adults2023, "c624f2e4228eea7296b2a77e2d4b7e50", "Env_exposure")
# Tyzzerella
source("Plot_specific_ASV_Abundance_Function.R")
plot.specific.ASV(adults2023, "80626c0d45293428d118ce1f05a1ab18", "Env_exposure")
