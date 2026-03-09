### Indicator Species Analysis ###

# Loading required libraries
library(tidyverse)
library(phyloseq)
library(indicspecies)

# Loading data
load("hiv_final.RData")

### Part 1) Indicator Species Analysis using ALL samples ###

# Creating IL-6 bins (High IL-6: >5 pg/mL & Low IL-6: ≤5 pg/mL)
IL6_binning_update <- data.frame(sample_data(hiv_final)) %>%
  mutate(IL6_bin = ifelse(IL.6_pg_mL > 5, "high", "low"))

# Replacing the existing metadata in the phyloseq object with the updated table that now includes IL-6 bin
sample_data(hiv_final) <- sample_data(IL6_binning_update)

# Determining the finest reliable taxonomic level available for aggregation
head(tax_table(hiv_final))
# Result: ASVs are resolved to the genus level, but not the species level

# Aggregating ASVs to the genus level
hiv_genus <- tax_glom(hiv_final, "Genus", NArm = FALSE)

# Converting counts to relative abundance
hiv_genus_RA <- transform_sample_counts(hiv_genus, fun = function(x) x / sum(x))

# Running the indicator species analysis
isa_il6 <- multipatt(t(otu_table(hiv_genus_RA)), cluster = sample_data(hiv_genus_RA)$IL6_bin)

# Looking at the summary of indicator taxa
summary(isa_il6)

# Creating taxonomy table
taxtable <- tax_table(hiv_genus_RA) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")

# Extracting significant indicator taxa
indicator_il6 <- isa_il6$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable) %>%
  filter(p.value <= 0.05)

View(indicator_il6)

# Creating a cleaner results table & labeling IL-6 group clearly
indicator_clean <- indicator_il6 %>%
  mutate(IL6_group = ifelse(s.high == 1, "High IL-6", "Low IL-6")) %>%
  select(Genus, Phylum, Family, stat, p.value, IL6_group) %>%
  arrange(desc(stat))

View(indicator_clean)

# Saving the results table
write.csv(indicator_clean, "indicator_species_IL6_all.csv", row.names = FALSE)


### Part 2) Indicator Species Analysis within PLWH only ###

# Subsetting dataset to PLWH (HIV+ samples only)
hiv_plwh <- subset_samples(hiv_final, HIV_Status == "Positive")

# Confirming the number of retained samples
nsamples(hiv_plwh)

# Removing taxa with zero counts after subsetting
hiv_plwh <- prune_taxa(taxa_sums(hiv_plwh) > 0, hiv_plwh)

# Checking IL-6 bin distribution with PLWH
table(sample_data(hiv_plwh)$IL6_bin)

# Aggregating taxa to genus level
hiv_plwh_genus <- tax_glom(hiv_plwh, "Genus", NArm = FALSE)

# Converting counts to relative abundance
hiv_plwh_genus_RA <- transform_sample_counts(hiv_plwh_genus, fun = function(x) x / sum(x))

# Running the indicator species analysis within PLWH
isa_il6_plwh <- multipatt(t(otu_table(hiv_plwh_genus_RA)), cluster = sample_data(hiv_plwh_genus_RA)$IL6_bin)

summary(isa_il6_plwh)

# Preparing taxonomy info for merging with ISA output
taxtable_plwh <- tax_table(hiv_plwh_genus_RA) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")

# Extracting significant taxa
indicator_il6_plwh <- isa_il6_plwh$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable_plwh) %>%
  filter(p.value <= 0.05)

View(indicator_il6_plwh)

# Creating a cleaner results table
indicator_clean_plwh <- indicator_il6_plwh %>%
  mutate(IL6_group = ifelse(s.high == 1, "High IL-6", "Low IL-6")) %>%
  select(Genus, Phylum, Family, stat, p.value, IL6_group) %>%
  arrange(desc(stat))

View(indicator_clean_plwh)

# Saving the PLWH results
write.csv(indicator_clean_plwh, "indicator_species_IL6_PLWH.csv", row.names = FALSE)




