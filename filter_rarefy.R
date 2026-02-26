# This is just for Elena's computer that is stuck in Chinese and I cannot read the console
# Sys.setenv(LANG="en")

library(phyloseq)
library(tidyverse)
library(vegan)

#### Filter #####
# Remove non-bacterial sequences, if any
hiv_filt <- subset_taxa(hiv,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
hiv_filt_nolow <- filter_taxa(hiv_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
hiv_filt_nolow_samps <- prune_samples(sample_sums(hiv_filt_nolow)>100, hiv_filt_nolow)
# Remove samples from Visit 3
hiv_final <- subset_samples(hiv_filt_nolow_samps, Visit_Cat == "2nd Visit")

#### Rarefy samples ####
# Generate rarefaction curve
# Saved as .png to repository
rarecurve(t(as.data.frame(otu_table(hiv_final))), cex=0.1)
# Used rngseed = 2 because we are Group 2 :P
# Sampling depth set to 25272 as per QIIME2 visualization
hiv_rare <- rarefy_even_depth(hiv_final, rngseed = 2, sample.size = 25272)
# 117 samples removed, 821 ASVs removed

#### Saving ####
save(hiv_final, file="hiv_final.RData")
save(hiv_rare, file="hiv_rare.RData")

# #### Other filters ####
# # Removed everything except for HIV status and IL-6 concentrations
# # Need to run "phyloseq object" Rscript first
# 
# samp_subset <- as.data.frame(meta[,c("HIV_Status", "Visit_Cat", "IL-6_pg_mL")]) # can edit to include only the variables we are interested in
# # Make sampleids the rownames
# rownames(samp_subset)<- meta$'sample-id'
# # Make phyloseq sample data with sample_data() function
# SAMP_sub <- sample_data(samp_subset)
# 
# hiv_subset <- phyloseq(OTU, SAMP_sub, TAX, phylotree)
# 
# # Remove non-bacterial sequences, if any
# hiv_subset_filt <- subset_taxa(hiv_subset,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# # Remove ASVs that have less than 5 counts total
# hiv_subset_filt_nolow <- filter_taxa(hiv_subset_filt, function(x) sum(x)>5, prune = TRUE)
# # Remove samples with less than 100 reads
# hiv_subset_filt_nolow_samps <- prune_samples(sample_sums(hiv_subset_filt_nolow)>100, hiv_subset_filt_nolow)
# # Remove samples from Visit 3
# hiv_subset_final <- subset_samples(hiv_subset_filt_nolow_samps, Visit_Cat == "2nd Visit")
# 
# # Rarefy
# hiv_subset_rare <- rarefy_even_depth(hiv_subset_final, rngseed = 2, sample.size = 25272)
# 
# save(hiv_subset_final, file="hiv_subset_final.RData")
# save(hiv_subset_rare, file="hiv_subset_rare.RData")