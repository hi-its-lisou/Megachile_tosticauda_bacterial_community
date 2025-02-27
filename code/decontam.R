################################################################################
#####                               DECONTAM                                ####
################################################################################
setwd("D:/Research/PhD/Manuscripts/Chapter 2/code/R project/Chapter2_analysis/data")

#load libraries
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

ps <- qza_to_phyloseq(
  features = "merged_table.qza",
  taxonomy = "merged_taxonomy.qza",
  tree = "merged_sepp_tree.qza",
  metadata = "combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

#DECONTAM

#Create a new column indicating which samples are negative controls
sample_data(ps)$is.neg <- sample_data(ps)$Negative_control == "TRUE"

#Run decontam prevalence method without threshold
ps_decontam_prev <- isContaminant(ps, method="prevalence", neg="is.neg")

##Now, let's plot the decontam scores to find a good threshold
ggplot(ps_decontam_prev, aes(x = p)) +
  geom_histogram(binwidth = 0.01) +
  labs(x = 'decontam Score', y='Number of species')

#Run decontam prevalence method with a 0.5 threshold
ps_decontam_prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)

#How many ASVs are considered contaminants with these thresholds?
table(ps_decontam_prev$contaminant)
table(ps_decontam_prev05$contaminant)

### Create prevalence/prevalence plots
#Create a new phyloseq object with table counts converted to presence/absence
ps_presence_absence <- transform_sample_counts(ps, function(abund) 1*(abund>0))

#Split phyloseq object into separate objects based on negative control vs. biological sample
ps_presence_absence_controls <- prune_samples(sample_data(ps)$Negative_control == "TRUE", ps_presence_absence)
ps_presence_absence_samples <- prune_samples(sample_data(ps)$Negative_control == "FALSE", ps_presence_absence)

#Merge them into a single dataframe, and add whether each ASV is considered a contaminant (default threshold)
ps_dataframe <- data.frame(samples = taxa_sums(ps_presence_absence_samples), 
                           controls = taxa_sums(ps_presence_absence_controls), 
                           contaminant = ps_decontam_prev$contaminant)

#Merge them into a single dataframe, and add whether each ASV is considered a contaminant (0.5 threshold)
ps_dataframe_05 <- data.frame(samples = taxa_sums(ps_presence_absence_samples), 
                              controls = taxa_sums(ps_presence_absence_controls), 
                              contaminant = ps_decontam_prev05$contaminant)
#Default threshold plot
ggplot(ps_dataframe, aes(x=controls, y=samples, color=contaminant)) + 
  geom_point(size=2) + 
  geom_jitter(height = 0.2, width = 0.2, size = 2) + 
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (Samples)") +
  guides(colour = guide_legend(title="Contaminant if >50% prevalence in negative controls", 
                               override.aes = list(size=10))) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "top",
  )

#0.5 threshold plot
ggplot(ps_dataframe_05, aes(x=controls, y=samples, color=contaminant)) + 
  geom_point(size=2) + 
  geom_jitter(height = 0.2, width = 0.2, size = 2) + 
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (Samples)") +
  guides(colour = guide_legend(title="Contaminant if >50% prevalence in negative controls", 
                               override.aes = list(size=10))) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "top",
  )

#Get taxonomy of putative contaminants
tax <- data.frame(as(tax_table(ps), "matrix"))

#Merge dataframes to obtain taxonomic information
df_merged <- merge(ps_decontam_prev05, tax, by=0, sort=FALSE)

#Generate two subsets of taxa: contaminants and non-contaminants
contaminant_ASVs <- subset(df_merged, contaminant == "TRUE")
non_contaminant_ASVs <- subset(df_merged, contaminant == "FALSE")

write_tsv(contaminant_ASVs, "putative_contaminants05.txt")
