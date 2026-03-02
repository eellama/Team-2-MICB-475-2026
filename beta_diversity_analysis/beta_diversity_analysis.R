#### Load libraries ####
library(phyloseq)
library(tidyverse)
library(vegan)

#### Load phyloseq object (RData) ####
load("filtered_and_rarefied/hiv_rare.RData")

# Make sure object is a phyloseq object
class(hiv_rare)

# Add column to indicate whether IL-6 levels are high or low into metadata table within phyloseq object
IL6_binning_update <- data.frame(sample_data(hiv_rare)) %>%
  mutate(IL6_bin = ifelse(IL.6_pg_mL > 5, "high", "low"))

sample_data(hiv_rare) <- sample_data(IL6_binning_update)

#### Beta diversity #####
# Create distance object for PCoA plot
wu_dist <- distance(hiv_rare, method="wunifrac")

# Create Weighted UniFrac PCoA plot
pcoa_wu <- ordinate(hiv_rare, method="PCoA", distance=wu_dist)

gg_pcoa <- plot_ordination(hiv_rare, pcoa_wu, color = "HIV_Status", shape="IL6_bin") +
  labs(pch="IL-6 levels (high/low)", col = "HIV status")
gg_pcoa

# Save plot
ggsave(filename = "beta_diversity_analysis/plot_pcoa.png",
       gg_pcoa,
       height=4, width=5)

### PERMANOVA (Permutational ANOVA) ####
# Calculate Weighted Unifrac distance matrix
dm_unifrac <- UniFrac(hiv_rare, weighted=TRUE)

# Plot dm_unifrac as an ordination to a PCoA plot
ord.unifrac <- ordinate(hiv_rare, method="PCoA", distance="unifrac")
plot_ordination(hiv_rare, ord.unifrac, color="HIV_Status", shape = "IL6_bin")

# Run the PERMANOVA on the above matrix for Weighted UniFrac
samp_dat_wdiv <- data.frame(sample_data(hiv_rare), estimate_richness(hiv_rare))
set.seed(500) # set.seed function is to ensure reproducibility of PERMANOVA results
adonis2(dm_unifrac ~ `HIV_Status`*IL6_bin, data=samp_dat_wdiv, by="terms")

# Confirm if significant difference in IL-6 bin is due to dispersion
bd <- betadisper(dm_unifrac, samp_dat_wdiv$IL6_bin)
anova(bd)

# Re-plot the original PCoA with ellipses to show a significant difference between IL-6 bin
gg_pcoa_ellipse <- plot_ordination(hiv_rare, ord.unifrac, color = "IL6_bin") +
  labs(color = "IL-6 bin") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.35,
           y = 0.3,
           label = "PERMANOVA\nR² = 0.027\np = 0.006")
gg_pcoa_ellipse

# Save plot
ggsave(filename = "beta_diversity_analysis/plot_pcoa_ellipse.png",
       gg_pcoa_ellipse,
       height=4, width=5)
