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
# Create distance matrix for PCoA plot
wu_dist <- distance(hiv_rare, method="wunifrac") # Weighted UniFrac
uu_dist <- distance(hiv_rare, method="unifrac") # Unweighted UniFrac
bc_dist <- distance(hiv_rare, method="bray") # Bray-Curtis
j_dist <- distance(hiv_rare, method = "jaccard", binary = TRUE) # Jaccard

# Create PCoA plot
pcoa_wu <- ordinate(hiv_rare, method="PCoA", distance=wu_dist)
pcoa_uu <- ordinate(hiv_rare, method="PCoA", distance=uu_dist)
pcoa_bc <- ordinate(hiv_rare, method="PCoA", distance=bc_dist)
pcoa_j <- ordinate(hiv_rare, method="PCoA", distance=j_dist)

gg_pcoa_wu <- plot_ordination(hiv_rare, pcoa_wu, color = "HIV_Status", shape="IL6_bin") +
  labs(pch="IL-6 levels (high/low)", col = "HIV status")
gg_pcoa_wu

gg_pcoa_uu <- plot_ordination(hiv_rare, pcoa_uu, color = "HIV_Status", shape="IL6_bin") +
  labs(pch="IL-6 levels (high/low)", col = "HIV status")
gg_pcoa_uu

gg_pcoa_bc <- plot_ordination(hiv_rare, pcoa_bc, color = "HIV_Status", shape="IL6_bin") +
  labs(pch="IL-6 levels (high/low)", col = "HIV status")
gg_pcoa_bc

gg_pcoa_j <- plot_ordination(hiv_rare, pcoa_j, color = "HIV_Status", shape="IL6_bin") +
  labs(pch="IL-6 levels (high/low)", col = "HIV status")
gg_pcoa_j

# Save plot
ggsave(filename = "beta_diversity_analysis/plots_all/plot_pcoa_wu.png",
       gg_pcoa_wu,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_all/plot_pcoa_uu.png",
       gg_pcoa_uu,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_all/plot_pcoa_bc.png",
       gg_pcoa_bc,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_all/plot_pcoa_j.png",
       gg_pcoa_j,
       height=4, width=5)

### PERMANOVA (Permutational ANOVA) ####
# Calculate distance matrix (dm)
dm_wu <- UniFrac(hiv_rare, weighted=TRUE) # Weighted UniFrac
dm_uu <- UniFrac(hiv_rare, weighted=FALSE) # Unweighted UniFrac
dm_bc <- vegdist(t(otu_table(hiv_rare)), method="bray") # Bray-Curtis
dm_j <- vegdist(t(otu_table(hiv_rare)), method="jaccard") # Jaccard

# Plot dm as an ordination to a PCoA plot
ord.wu <- ordinate(hiv_rare, method="PCoA", distance="wunifrac")
plot_ordination(hiv_rare, ord.wu, color="HIV_Status", shape = "IL6_bin")

ord.uu <- ordinate(hiv_rare, method="PCoA", distance="unifrac")
plot_ordination(hiv_rare, ord.uu, color="HIV_Status", shape = "IL6_bin")

ord.bc <- ordinate(hiv_rare, method="PCoA", distance="bray")
plot_ordination(hiv_rare, ord.bc, color="HIV_Status", shape = "IL6_bin")

ord.j <- ordinate(hiv_rare, method="PCoA", distance="jaccard")
plot_ordination(hiv_rare, ord.j, color="HIV_Status", shape = "IL6_bin")

# Run the PERMANOVA on the above matrix for Weighted UniFrac
samp_dat_wdiv <- data.frame(sample_data(hiv_rare), estimate_richness(hiv_rare))
set.seed(500) # set.seed function is to ensure reproducibility of PERMANOVA results
adonis2(dm_wu ~ `HIV_Status`*IL6_bin, data=samp_dat_wdiv, by="terms")
adonis2(dm_uu ~ `HIV_Status`*IL6_bin, data=samp_dat_wdiv, by="terms")
adonis2(dm_bc ~ `HIV_Status`*IL6_bin, data=samp_dat_wdiv, by="terms")
adonis2(dm_j ~ `HIV_Status`*IL6_bin, data=samp_dat_wdiv, by="terms")

# Confirm if significant difference in IL-6 bin is due to dispersion
bd_wu <- betadisper(dm_wu, samp_dat_wdiv$IL6_bin)
anova(bd_wu)

bd_uu <- betadisper(dm_uu, samp_dat_wdiv$HIV_Status)
anova(bd_uu)

bd_bc <- betadisper(dm_bc, samp_dat_wdiv$HIV_Status)
anova(bd_bc)

bd_j <- betadisper(dm_j, samp_dat_wdiv$HIV_Status)
anova(bd_j)

# Re-plot the original PCoA with ellipses to show a significant difference between IL-6 bin or HIV status
gg_pcoa_wu_ellipse <- plot_ordination(hiv_rare, ord.wu, color = "IL6_bin") +
  labs(color = "IL-6 bin") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.050,
           y = -0.03,
           label = "PERMANOVA\nR² = 0.027\np = 0.006")
gg_pcoa_wu_ellipse

gg_pcoa_uu_ellipse <- plot_ordination(hiv_rare, ord.uu, color = "HIV_Status") +
  labs(color = "HIV status") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.40,
           y = 0.25,
           label = "PERMANOVA\nR² = 0.016\np = 0.017")
gg_pcoa_uu_ellipse

gg_pcoa_bc_ellipse <- plot_ordination(hiv_rare, ord.bc, color = "HIV_Status") +
  labs(color = "HIV status") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.4,
           y = 0.3,
           label = "PERMANOVA\nR² = 0.018\np = 0.020")
gg_pcoa_bc_ellipse

gg_pcoa_j_ellipse <- plot_ordination(hiv_rare, ord.j, color = "HIV_Status") +
  labs(color = "HIV status") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.35,
           y = 0.3,
           label = "PERMANOVA\nR² = 0.016\np = 0.021")
gg_pcoa_j_ellipse

# Save plot
ggsave(filename = "beta_diversity_analysis/plots_all/plot_pcoa_wu_ellipse.png",
       gg_pcoa_wu_ellipse,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_all/plot_pcoa_uu_ellipse.png",
       gg_pcoa_uu_ellipse,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_all/plot_pcoa_bc_ellipse.png",
       gg_pcoa_bc_ellipse,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_all/plot_pcoa_j_ellipse.png",
       gg_pcoa_j_ellipse,
       height=4, width=5)

#### Identifying which participants the 3 diverging dots identified in wu plots belong to ####
# Extract coordinates
pcoa_df <- plot_ordination(hiv_rare, ord.wu, justDF = TRUE)

# Subset participants of interest
pcoa_df_subset <- subset(pcoa_df, Axis.1 > 0.025 & Axis.2 > 0.01) %>%
  select(Axis.1, Axis.2, HIV_Status, IL.6_pg_mL, IL6_bin)
view(pcoa_df_subset)

# Convert metadata in phyloseq object hiv_rare to df
metadata <- data.frame(sample_data(hiv_rare)) %>%
  rownames_to_column(var = "SampleID")

nrow(metadata)

# See how many samples excluding the 3 identified samples have IL-6 above 17 pg/mL
metadata_exclude_3 <- metadata %>%
  filter(!SampleID %in% c("ERR12057704", "ERR12063224", "ERR12063266") &
           IL.6_pg_mL >= 17) %>%
  select(SampleID, HIV_Status, IL.6_pg_mL, IL6_bin)

nrow(metadata_exclude_3) # Highest IL-6 lvls

# All samples with IL-6 > 17 pg/mL
metadata_above17 <- metadata %>%
  filter(IL.6_pg_mL > 17) %>%
  select("SampleID", "HIV_Status", "IL.6_pg_mL", "IL6_bin") %>%
  arrange(desc(IL.6_pg_mL))

# Summary report
## 87 participants total (both HIV+ and HIV-)
## 3 participants (2 HIV+, 1 HIV-) identified from PCoA plot. All have IL-6 levels > 17 pg/mL
## Excluding the above 3, 6 participants have IL-6 levels > 10 pg/mL. All HIV+