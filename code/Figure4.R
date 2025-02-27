# Load the required libraries ####
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
library(patchwork)

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

# Subset the data for samples that underwent qPCRs (have CT scores)
qPCRs <- subset_samples(ps, complete.cases(AVG_CTscore) &
                          sample_type %in% c("Negative_control", "Food", "Adults", "Honey_bee"))

# Make the phyloseq into a dataframe
sample_data_df <- as.data.frame(sample_data(qPCRs))

# Custom order for figure
custom_order <- c("Negative_control", "Adults", "Food", "Honey_bee")
sample_data_df$sample_type <- factor(sample_data_df$sample_type, levels = custom_order)

# Convert NA values to strings so that they are still included in the plot
sample_data_df$Env_exposure <- as.character(sample_data_df$Env_exposure)
sample_data_df$Env_exposure[is.na(sample_data_df$Env_exposure)] <- "NA"
sample_data_df$Env_exposure <- factor(sample_data_df$Env_exposure, levels = c("Free-flying", "Nest_only", "None", "NA"))

# Custom colors for sample types
custom_colors <- c("Negative_control" = "red", "Adults" = "black", "Food" = "yellow", "Honey_bee" = "blue")

# Create the ggplot
qPCRPlot <- ggplot(sample_data_df, aes(x = sampleid, y = logDNA, color = sample_type)) +
  geom_hline(aes(linetype = "limit of detection", yintercept = 2.339477), color = "green") +
  geom_point(aes(shape = Env_exposure), size = 3, na.rm = TRUE) +
  scale_shape_manual(values = c(
    "Free-flying" = 1,
    "Nest_only" = 4,
    "None" = 3,
    "NA" = 16
  )) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(name = "logDNA") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size = 12)) +
  facet_wrap(~sample_type, ncol = 5, scale = "free_x")

qPCRPlot

# Save the figure 
ggsave("figures/OctFigure4.png", height=5, width=10)
   
   