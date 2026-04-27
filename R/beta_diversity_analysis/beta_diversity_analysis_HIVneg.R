#### Load libraries ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(showtext)

# Add Calibri font
font_add(family = "Calibri", regular = "Calibri.ttf")
showtext_auto()

#### Load phyloseq object (RData) ####
load("filtered_and_rarefied/hiv_rare.RData")

# Make sure object is a phyloseq object
class(hiv_rare)

# Add column to indicate whether IL-6 levels are high or low into metadata table within phyloseq object
IL6_binning_update <- data.frame(sample_data(hiv_rare)) %>%
  mutate(IL6_bin = ifelse(IL.6_pg_mL > 5, "high", "low"))

sample_data(hiv_rare) <- sample_data(IL6_binning_update)

# Subset HIV- participants only
hivneg_only_rare <- subset_samples(hiv_rare, HIV_Status == "Negative")

#### Beta diversity #####
# Create distance matrix for PCoA plot
wu_dist <- distance(hivneg_only_rare, method="wunifrac") # Weighted UniFrac
uu_dist <- distance(hivneg_only_rare, method="unifrac") # Unweighted UniFrac
bc_dist <- distance(hivneg_only_rare, method="bray") # Bray-Curtis
j_dist <- distance(hivneg_only_rare, method = "jaccard", binary = TRUE) # Jaccard

# Create PCoA plot
pcoa_wu <- ordinate(hivneg_only_rare, method="PCoA", distance=wu_dist)
pcoa_uu <- ordinate(hivneg_only_rare, method="PCoA", distance=uu_dist)
pcoa_bc <- ordinate(hivneg_only_rare, method="PCoA", distance=bc_dist)
pcoa_j <- ordinate(hivneg_only_rare, method="PCoA", distance=j_dist)

gg_pcoa_wu <- plot_ordination(hivneg_only_rare, pcoa_wu, color = "IL6_bin") +
  labs(col="IL-6 levels")
gg_pcoa_wu

gg_pcoa_uu <- plot_ordination(hivneg_only_rare, pcoa_uu, color = "IL6_bin") +
  labs(col="IL-6 levels")
gg_pcoa_uu

gg_pcoa_bc <- plot_ordination(hivneg_only_rare, pcoa_bc, color = "IL6_bin") +
  labs(col="IL-6 levels")
gg_pcoa_bc

gg_pcoa_j <- plot_ordination(hivneg_only_rare, pcoa_j, color = "IL6_bin") +
  labs(col="IL-6 levels")
gg_pcoa_j

# Save plot
ggsave(filename = "beta_diversity_analysis/plots_hivneg/plot_pcoa_wu_HIVneg.png",
       gg_pcoa_wu,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_hivneg/plot_pcoa_uu_HIVneg.png",
       gg_pcoa_uu,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_hivneg/plot_pcoa_bc_HIVneg.png",
       gg_pcoa_bc,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_hivneg/plot_pcoa_j_HIVneg.png",
       gg_pcoa_j,
       height=4, width=5)

### PERMANOVA (Permutational ANOVA) ####
# Calculate distance matrix (dm)
dm_wu <- UniFrac(hivneg_only_rare, weighted=TRUE) # Weighted UniFrac
dm_uu <- UniFrac(hivneg_only_rare, weighted=FALSE) # Unweighted UniFrac
dm_bc <- vegdist(t(otu_table(hivneg_only_rare)), method="bray") # Bray-Curtis
dm_j <- vegdist(t(otu_table(hivneg_only_rare)), method="jaccard") # Jaccard

# Plot dm as an ordination to a PCoA plot
ord.wu <- ordinate(hivneg_only_rare, method="PCoA", distance="wunifrac")
plot_ordination(hivneg_only_rare, ord.wu, color="IL6_bin")

ord.uu <- ordinate(hivneg_only_rare, method="PCoA", distance="unifrac")
plot_ordination(hivneg_only_rare, ord.uu, color="IL6_bin")

ord.bc <- ordinate(hivneg_only_rare, method="PCoA", distance="bray")
plot_ordination(hivneg_only_rare, ord.bc, color="IL6_bin")

ord.j <- ordinate(hivneg_only_rare, method="PCoA", distance="jaccard")
plot_ordination(hivneg_only_rare, ord.j, color="IL6_bin")

# Run the PERMANOVA on the above matrix## Only IL6_bin as response variable?
samp_dat_wdiv <- data.frame(sample_data(hivneg_only_rare), estimate_richness(hivneg_only_rare))
set.seed(500) # set.seed function is to ensure reproducibility of PERMANOVA results
adonis2(dm_wu ~ IL6_bin, data=samp_dat_wdiv, by="terms")
adonis2(dm_uu ~ IL6_bin, data=samp_dat_wdiv, by="terms")
adonis2(dm_bc ~ IL6_bin, data=samp_dat_wdiv, by="terms")
adonis2(dm_j ~ IL6_bin, data=samp_dat_wdiv, by="terms")

# Re-plot the original PCoA with ellipses
gg_pcoa_wu_ellipse <- plot_ordination(hivneg_only_rare, ord.wu, color = "IL6_bin") +
  labs(title = "HIV- participants", color = "IL-6 bin", x = "PC1(20%)", y = "PC2(13.8%)") +
  scale_color_manual(values = c(
    "low" = "#1f78b4",
    "high" = "#e31a1c")) +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.045,
           y = 0.025,
           label = "PERMANOVA\nR² = 0.024\np = 0.571",
           size = 5,
           lineheight = 0.8,
           family = "Calibri") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(family = "Calibri", size = 17))
gg_pcoa_wu_ellipse

gg_pcoa_uu_ellipse <- plot_ordination(hivneg_only_rare, ord.uu, color = "IL6_bin") +
  labs(color = "IL-6 bin") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.40,
           y = 0.25,
           label = "PERMANOVA\nR² = 0.032\np = 0.332")
gg_pcoa_uu_ellipse

gg_pcoa_bc_ellipse <- plot_ordination(hivneg_only_rare, ord.bc, color = "IL6_bin") +
  labs(color = "IL-6 bin") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.4,
           y = 0.3,
           label = "PERMANOVA\nR² = 0.029\np = 0.554")
gg_pcoa_bc_ellipse

gg_pcoa_j_ellipse <- plot_ordination(hivneg_only_rare, ord.j, color = "IL6_bin") +
  labs(color = "IL-6 bin") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = -0.45,
           y = 0.3,
           label = "PERMANOVA\nR² = 0.030\np = 0.524")
gg_pcoa_j_ellipse

# Save plot
ggsave(filename = "beta_diversity_analysis/plots_hivneg/plot_pcoa_wu_ellipse_HIVneg.png",
       gg_pcoa_wu_ellipse,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_hivneg/plot_pcoa_uu_ellipse_HIVneg.png",
       gg_pcoa_uu_ellipse,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_hivneg/plot_pcoa_bc_ellipse_HIVneg.png",
       gg_pcoa_bc_ellipse,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_hivneg/plot_pcoa_j_ellipse_HIVneg.png",
       gg_pcoa_j_ellipse,
       height=4, width=5)
