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
library(scales)

# Load in sequencing data
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

################################################################################
#####          calculating the % of chloroplasts and mitochondria           ####
################################################################################

# Function to computer percentage
compute_percentage <- function(taxon, ps) {
  counts_data <- as.data.frame(otu_table(ps))  #Extract the counts data from the phyloseq object
  taxon_row <- which(tax_table(ps)[, "Genus"] == taxon)  #Identify the row corresponding to the specific taxon in the taxonomic table
  taxon_counts <- sum(counts_data[taxon_row, ])  #Sum the counts of the specific taxon across all samples
  total_counts <- sum(counts_data)  #Calculate the percentage of reads attributed to the specific taxon
  percentage_taxon <- (taxon_counts / total_counts) * 100
  return(percentage_taxon)
}

# Sublist sample types
Prepupae <- subset_samples(ps, sample_type %in% c("Prepupae")& !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic")))
Larva <- subset_samples(ps, sample_type %in% c("Larvae"))
Honey_bee <- subset_samples(ps, sample_type %in% c("Honey_bee"))
Food <- subset_samples(ps, sample_type %in% c("Food"))
AdultsFree <- subset_samples(ps, Env_exposure %in% c("Free-flying")) #seperate only natural adult samples

AdultsNest <- subset_samples(ps, Env_exposure %in% c("Nest_only"))
AdultsNone <- subset_samples(ps, Env_exposure %in% c("None"))

# Run function for each sublisted sample type
# For Chloroplast
percentage_Chloroplast_Adults <- compute_percentage("Chloroplast", AdultsFree) 
percentage_Chloroplast_Larva <- compute_percentage("Chloroplast", Larva)
percentage_Chloroplast_Prepupae <- compute_percentage("Chloroplast", Prepupae)
percentage_Chloroplast_Food <- compute_percentage("Chloroplast", Food)
percentage_Chloroplast_HB <- compute_percentage("Chloroplast", Honey_bee)

percentage_Chloroplast_Adults
percentage_Chloroplast_Food
percentage_Chloroplast_Larva
percentage_Chloroplast_Prepupae
percentage_Chloroplast_HB

percentage_Chloroplast_AdultsNest <- compute_percentage("Chloroplast", AdultsNest)
percentage_Chloroplast_AdultsNone <- compute_percentage("Chloroplast", AdultsNone)

percentage_Chloroplast_antibiotic
percentage_Chloroplast_control
percentage_Chloroplast_AdultsNest
percentage_Chloroplast_AdultsNone

# For Mitochondria
percentage_Mitochondria_Adults <- compute_percentage("Mitochondria", AdultsFree) 
percentage_Mitochondria_Larva <- compute_percentage("Mitochondria", Larva)
percentage_Mitochondria_Prepupae <- compute_percentage("Mitochondria", Prepupae)
percentage_Mitochondria_Food <- compute_percentage("Mitochondria", Food)
percentage_Mitochondria_HB <- compute_percentage("Mitochondria", Honey_bee)

percentage_Mitochondria_Adults
percentage_Mitochondria_Food
percentage_Mitochondria_Larva
percentage_Mitochondria_Prepupae
percentage_Mitochondria_HB

percentage_Mitochondria_AdultsNest <- compute_percentage("Mitochondria", AdultsNest)
percentage_Mitochondria_AdultsNone <- compute_percentage("Mitochondria", AdultsNone)

percentage_Mitochondria_AdultsNest
percentage_Mitochondria_AdultsNone

### Calculating the proportion of putitive contaminate asvs (called using decontam at a threshold of 0.5)
#vector for contaminant asvs
contam_asvs <- read_delim("input_files/decontam_asvs.txt", 
                                        delim = "\n", 
                                        col_names = "asv")
contam_asvs <- contam_asvs$asv
all_asvs <- taxa_names(ps)
asvs_to_keep <- all_asvs[!(all_asvs %in% contam_asvs)]
ps_no_contam <- prune_taxa(asvs_to_keep, ps)
ps_no_contam #pyloseq object without any contaminant asvs

#Adult tostis contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Free-flying")))))) #Specifically for free-fying the adults
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, Env_exposure %in% c("Free-flying"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Pollen provision contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Food")))))) 
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Food")))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Larvae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Larvae"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Larvae")))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Prepupae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == "Prepupae" & !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic")))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Prepupae" & !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic"))))))) 
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Honey bee contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Honey_bee")))))) 
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Honey_bee"))))) 
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Nest eclosed bees contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Nest_only"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, Env_exposure %in% c("Nest_only"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Controlled eclosion bees contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("None"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, Env_exposure %in% c("None"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

### Working biomass after removing ALL non-target reads ###
# Remove contaminants, chloroplasts, and mitochondria ####
chloro_mito_decontam_asvs <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                     delim = "\n", 
                     col_names = "asv")
chloro_mito_decontam_asvs <- chloro_mito_decontam_asvs$asv
all_asvs <- taxa_names(ps)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps_filtered <- prune_taxa(asvs_to_keep, ps)
ps_filtered #phyloseq object with all off-target reads removed

# Free flying adults working biomass 
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Free-flying"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered,Env_exposure %in% c("Free-flying"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#Larvae working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Larvae")))))) #replace with sample_type in question
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Larvae"))))) #replace with sample_type in question
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#Prepupae working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Prepupae") & !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic"))))))) 
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Prepupae" & !(sample_data(ps_filtered)$AB_treatment %in% c("control", "antibiotic")))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#Food working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Food")))))) #replace with sample_type in question
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Food"))))) #replace with sample_type in question
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#honey bee working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Honey_bee")))))) #replace with sample_type in question
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Honey_bee"))))) #replace with sample_type in question
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

# nest eclosed adults working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Nest_only"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, Env_exposure %in% c("Nest_only"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

# controlled eclosion adults working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("None"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, Env_exposure %in% c("None"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

