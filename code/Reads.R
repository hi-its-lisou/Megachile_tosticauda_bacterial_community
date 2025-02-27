# Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(vegan)
library(ggplot2)
library(dplyr)


### Load raw, unfiltered data object ####
ps_raw <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_5.tsv")
ps_raw@sam_data$sampleid = rownames(ps_raw@sam_data)
ps_raw


### Use only the samples from manuscript
# Honey bees
HB_subset <- subset_samples(ps_raw, sample_type %in% c("Honey_bee"))
# Prepupae
prepupal_subset <- subset_samples(ps_raw, sample_type %in% c("Prepupae") & is.na(AB_treatment) & !Nest_id %in% c("5", "13", "14", "27"))
# Larvae
larva_subset <- subset_samples(ps_raw, sample_type %in% c("Larvae") & Year %in% c("2022"))
# Adults
adults <- subset_samples(ps_raw, sample_type %in% c("Adults"))
adults
# Frass
frass_subset <- subset_samples(ps_raw, Description_of_sample %in% c("CocoonFrass", "Frass", "PinkCocoonFrass"))
frass_subset@sam_data$sample_type[which(frass_subset@sam_data$sample_type == "Nest_contents")] <- "Frass contents"
# Pollen provisions
food_through_time <- subset_samples(ps_raw, sample_type %in% c("Food") 
                                    & Nest_id %in% c("5", "13", "14", "27", "15")
                                    & Cell_ID %in% c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A"))
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "5")] <- "1"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "13")] <- "2"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "14")] <- "3"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "27")] <- "4"
# Negative controls
negative_controls <- (subset_samples(ps_raw, sample_type %in% c("Negative_control")))

# Merge all together
ps_raw <- merge_phyloseq(food_through_time, larva_subset, prepupal_subset, adults, HB_subset, frass_subset, negative_controls)


### Total reads in the raw (unfiltered) dataset
total_reads_raw <- sum(sample_sums(ps_raw))
sort(sample_sums(ps_raw))


### Remove contaminants, chloroplasts, and mitochondria ####
contam <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                     delim = "\n", 
                     col_names = "asv")
chloro_mito_decontam_asvs <- contam$asv
all_asvs <- taxa_names(ps_raw)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps_filtered <- prune_taxa(asvs_to_keep, ps_raw)
ps_filtered

sort(sample_sums(ps_filtered))

# Total reads in the filtered dataset
total_reads_filtered <- sum(sample_sums(ps_filtered))

# Reads per sample in the raw dataset
reads_per_sample_raw <- sample_sums(ps_raw)
mean_reads_raw <- mean(reads_per_sample_raw)

# Reads per sample in the filtered dataset
reads_per_sample_filtered <- sample_sums(ps_filtered)
mean_reads_filtered <- mean(reads_per_sample_filtered)

# Summary
cat("Sequencing Summary:\n")
cat("Total reads generated (raw):", total_reads_raw, "\n")
cat("Total reads remaining after filtering:", total_reads_filtered, "\n")
cat("Mean reads per sample (raw):", round(mean_reads_raw), "\n")
cat("Mean reads per sample (filtered):", round(mean_reads_filtered), "\n")


### Making a rarefaction curve for filtered data ###
# Remove contaminants, chloroplasts, and mitochondria
mat <- t(otu_table(filtered_physeq))
raremax <- min(rowSums(mat))

taxa_are_rows(filtered_physeq)
mat <- t(otu_table(filtered_physeq))
class(mat) <- "matrix"

#  Setting class(x) to "matrix" sets attribute to NULL; result will no longer be an S4 object
class(mat)

mat <- as(t(otu_table(filtered_physeq)), "matrix")
class(mat)

raremax <- min(rowSums(mat))

# Plot the curve with a max sequencing depth of 2500
rarecurve(mat, step = 100, 
          sample = raremax, 
          col = "blue", 
          label = FALSE, 
          xlim = c(0, 2500),  
          xlab = "Sequencing Depth", 
          ylab = "Observed Species")

### Calculate the proportion of reads that were off-target ###
# Extract the OTU table and taxonomy table from the phyloseq object
otu_table_df <- as.data.frame(as(otu_table(ps_raw), "matrix"))
tax_table_df <- as.data.frame(as(tax_table(ps_raw), "matrix"))

# Ensure the taxa names are consistent in both tables
otu_table_df$TaxaID <- rownames(otu_table_df)
tax_table_df$TaxaID <- rownames(tax_table_df)

# Merge OTU and taxonomy tables based on TaxaID
merged_data <- merge(otu_table_df, tax_table_df, by = "TaxaID")

# Identify chloroplast and mitochondrial ASVs/OTUs in the taxonomy
chloroplast_mitochondria_taxa <- merged_data %>%
  filter(grepl("chloroplast", Genus, ignore.case = TRUE) | 
           grepl("mitochondria", Family, ignore.case = TRUE))

# Select only the OTU count columns (excluding non-numeric columns like taxonomic information)
otu_columns <- select(chloroplast_mitochondria_taxa, where(is.numeric))

# Calculate total reads for chloroplast and mitochondria ASVs/OTUs
chloroplast_mitochondria_reads <- rowSums(otu_columns)

# Calculate total reads in the entire dataset
total_reads <- sum(rowSums(otu_table(ps_raw)))

# Calculate the proportion of reads that are chloroplasts or mitochondria
proportion_chloroplast_mitochondria <- sum(chloroplast_mitochondria_reads) / total_reads

# Print the results
cat("Total reads:", total_reads, "\n")
cat("Chloroplast/Mitochondria reads:", sum(chloroplast_mitochondria_reads), "\n")
cat("Proportion of reads that are chloroplast or mitochondria:", round(proportion_chloroplast_mitochondria * 100, 2), "%\n")
<<<<<<< HEAD


mat <- t(otu_table(ps_raw))
raremax <- min(rowSums(mat))

#> Loading required package: permute
#> Loading required package: lattice
#> This is vegan 2.6-4

taxa_are_rows(ps_raw)
mat <- t(otu_table(ps_raw))
class(mat) <- "matrix"
class(mat)
#> [1] "matrix" "array"

mat <- as(t(otu_table(ps_raw)), "matrix")
#> [1] "matrix" "array"

raremax <- min(rowSums(mat))
rarecurve(mat, step = 20, sample = raremax, col = "blue", label = FALSE, cex = 0.6, xlab = "Sequencing depth", ylab = "Observed Species")
=======
>>>>>>> c90a874496f061303cb99cf4d6e4683ee5e0831f
