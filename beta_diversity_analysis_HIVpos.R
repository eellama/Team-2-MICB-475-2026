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

# Subset HIV+ participants only
hivpos_only_rare <- subset_samples(hiv_rare, HIV_Status == "Positive")

#### Beta diversity #####
# Create distance matrix for PCoA plot
wu_dist <- distance(hivpos_only_rare, method="wunifrac") # Weighted UniFrac
uu_dist <- distance(hivpos_only_rare, method="unifrac") # Unweighted UniFrac
bc_dist <- distance(hivpos_only_rare, method="bray") # Bray-Curtis
j_dist <- distance(hivpos_only_rare, method = "jaccard", binary = TRUE) # Jaccard

# Create PCoA plot
pcoa_wu <- ordinate(hivpos_only_rare, method="PCoA", distance=wu_dist)
pcoa_uu <- ordinate(hivpos_only_rare, method="PCoA", distance=uu_dist)
pcoa_bc <- ordinate(hivpos_only_rare, method="PCoA", distance=bc_dist)
pcoa_j <- ordinate(hivpos_only_rare, method="PCoA", distance=j_dist)

gg_pcoa_wu <- plot_ordination(hivpos_only_rare, pcoa_wu, color = "IL6_bin") +
  labs(col="IL-6 levels")
gg_pcoa_wu

gg_pcoa_uu <- plot_ordination(hivpos_only_rare, pcoa_uu, color = "IL6_bin") +
  labs(col="IL-6 levels")
gg_pcoa_uu

gg_pcoa_bc <- plot_ordination(hivpos_only_rare, pcoa_bc, color = "IL6_bin") +
  labs(col="IL-6 levels")
gg_pcoa_bc

gg_pcoa_j <- plot_ordination(hivpos_only_rare, pcoa_j, color = "IL6_bin") +
  labs(col="IL-6 levels")
gg_pcoa_j

# Save plot
ggsave(filename = "beta_diversity_analysis/plots_hivpos/plot_pcoa_wu_HIVpos.png",
       gg_pcoa_wu,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_hivpos/plot_pcoa_uu_HIVpos.png",
       gg_pcoa_uu,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_hivpos/plot_pcoa_bc_HIVpos.png",
       gg_pcoa_bc,
       height=4, width=5)

ggsave(filename = "beta_diversity_analysis/plots_hivpos/plot_pcoa_j_HIVpos.png",
       gg_pcoa_j,
       height=4, width=5)

### PERMANOVA (Permutational ANOVA) ####
# Calculate distance matrix (dm)
dm_wu <- UniFrac(hivpos_only_rare, weighted=TRUE) # Weighted UniFrac
dm_uu <- UniFrac(hivpos_only_rare, weighted=FALSE) # Unweighted UniFrac
dm_bc <- vegdist(t(otu_table(hivpos_only_rare)), method="bray") # Bray-Curtis
dm_j <- vegdist(t(otu_table(hivpos_only_rare)), method="jaccard") # Jaccard

# Plot dm as an ordination to a PCoA plot
## Plot just color="IL6_bin"?
ord.wu <- ordinate(hivpos_only_rare, method="PCoA", distance="wunifrac")
plot_ordination(hivpos_only_rare, ord.wu, color="IL6_bin")

ord.uu <- ordinate(hivpos_only_rare, method="PCoA", distance="unifrac")
plot_ordination(hivpos_only_rare, ord.uu, color="IL6_bin")

ord.bc <- ordinate(hivpos_only_rare, method="PCoA", distance="bray")
plot_ordination(hivpos_only_rare, ord.bc, color="IL6_bin")

ord.j <- ordinate(hivpos_only_rare, method="PCoA", distance="jaccard")
plot_ordination(hivpos_only_rare, ord.j, color="IL6_bin")

# Run the PERMANOVA on the above matrix
## Only IL6_bin as response variable?
samp_dat_wdiv <- data.frame(sample_data(hivpos_only_rare), estimate_richness(hivpos_only_rare))
set.seed(500) # set.seed function is to ensure reproducibility of PERMANOVA results
adonis2(dm_wu ~ IL6_bin, data=samp_dat_wdiv, by="terms")
adonis2(dm_uu ~ IL6_bin, data=samp_dat_wdiv, by="terms") # p-value insignificant
adonis2(dm_bc ~ IL6_bin, data=samp_dat_wdiv, by="terms") # p-value insignificant
adonis2(dm_j ~ IL6_bin, data=samp_dat_wdiv, by="terms") # p-value insignificant

# Confirm if significant difference in IL-6 bin is due to dispersion
bd_wu <- betadisper(dm_wu, samp_dat_wdiv$IL6_bin)
anova(bd_wu)

# Re-plot the original PCoA with ellipses
gg_pcoa_wu_ellipse <- plot_ordination(hivpos_only_rare, ord.wu, color = "IL6_bin") +
  labs(color = "IL-6 bin") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.045,
           y = -0.025,
           label = "PERMANOVA\nR² = 0.041\np = 0.01")
gg_pcoa_wu_ellipse

gg_pcoa_uu_ellipse <- plot_ordination(hivpos_only_rare, ord.uu, color = "IL6_bin") +
  labs(color = "IL-6 bin") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.37,
           y = -0.3,
           label = "PERMANOVA\nR² = 0.019\np = 0.436")
gg_pcoa_uu_ellipse

gg_pcoa_bc_ellipse <- plot_ordination(hivpos_only_rare, ord.bc, color = "IL6_bin") +
  labs(color = "IL-6 bin") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.4,
           y = 0.3,
           label = "PERMANOVA\nR² = 0.022\np = 0.244")
gg_pcoa_bc_ellipse

gg_pcoa_j_ellipse <- plot_ordination(hivpos_only_rare, ord.j, color = "IL6_bin") +
  labs(color = "IL-6 bin") +
  stat_ellipse(type = "norm") +
  annotate("text",
           x = 0.42,
           y = 0.3,
           label = "PERMANOVA\nR² = 0.020\np = 0.251")
gg_pcoa_j_ellipse

# Save plot
ggsave(filename = "beta_diversity_analysis/plots_hivpos/plot_pcoa_wu_ellipse_HIVpos.png",
       gg_pcoa_wu_ellipse,
       height=4, width=6)

ggsave(filename = "beta_diversity_analysis/plots_hivpos/plot_pcoa_uu_ellipse_HIVpos.png",
       gg_pcoa_uu_ellipse,
       height=4, width=6)

ggsave(filename = "beta_diversity_analysis/plots_hivpos/plot_pcoa_bc_ellipse_HIVpos.png",
       gg_pcoa_bc_ellipse,
       height=4, width=6)

ggsave(filename = "beta_diversity_analysis/plots_hivpos/plot_pcoa_j_ellipse_HIVpos.png",
       gg_pcoa_j_ellipse,
       height=4, width=6)

#### Identifying which participants the 3 diverging dots identified in wu plots belong to ####
# Extract coordinates
pcoa_df <- plot_ordination(hivpos_only_rare, ord.wu, justDF = TRUE)

# Subset participants of interest
pcoa_df_subset <- subset(pcoa_df, Axis.1 < -0.03 & Axis.2 < 0) %>%
  select(Axis.1, Axis.2, HIV_Status, IL.6_pg_mL, IL6_bin)
view(pcoa_df_subset)

# Extract sample names
clustered_samples <- rownames(pcoa_df_subset)

# Summary report
## 2 participants identified from PCoA plot
## Both are present among the 3 participants identified in HIVall

# Create dot plot of IL-6 lvl distribution with the 2 identified samples highlighted
dist_meta <- data.frame(sample_data(hivpos_only_rare)) %>%
  select("IL.6_pg_mL", "IL6_bin")

dist_meta_abv17 <- dist_meta %>%
  rownames_to_column("sample_name") %>%
  filter(IL.6_pg_mL > 17) %>% 
  arrange(desc(IL.6_pg_mL))

# Prepare dist_meta for plotting
dist_meta$highlight <- ifelse(rownames(dist_meta) %in% clustered_samples,
                         "Outlier",
                         "Normal")

dist_plot_highlight <- ggplot(dist_meta, aes(x = IL.6_pg_mL, fill = highlight)) +
  geom_dotplot(binwidth = 0.2, stackdir = "up", dotsize = 1.2, stackratio = 1.5) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  scale_fill_manual(name = "PCoA status", values = c("grey70", "red")) +
  labs(x = "IL-6 levels (pg/mL)", y = NULL) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 0.2)

# Save plot
ggsave(filename = "beta_diversity_analysis/IL6_distplot_highlight.png",
       dist_plot_highlight)

# Check where other samples with IL-6 levels > 17 pg/mL are on PCoA plot
coords <- as.data.frame(ord.wu$vectors) %>%
  rownames_to_column("sample_name") %>%
  select(sample_name, Axis.1, Axis.2)
  
coords_filt <- coords[coords$sample_name %in% dist_meta_abv17$sample_name, ]

dist_meta_abv17_coords <- left_join(dist_meta_abv17, coords_filt, by = "sample_name") %>%
  arrange(desc(IL.6_pg_mL))

# Recreate gg_pcoa_wu_ellipse with arrow pointing at sample with highest IL-6 levels
highest_IL6 <- dist_meta_abv17_coords[1, ]

gg_pcoa_wu_ellipse_arrow <- gg_pcoa_wu_ellipse +
  annotate(
    "segment",
    x = highest_IL6$Axis.1 - 0.005,
    y = highest_IL6$Axis.2 + 0.005,
    xend = highest_IL6$Axis.1 - 0.001,
    yend = highest_IL6$Axis.2 + 0.001,
    arrow = arrow(length = unit(0.2, "cm"))
  )

# Save plot
ggsave(filename = "beta_diversity_analysis/plots_hivpos/plot_pcoa_wu_ellipse_arrow_HIVpos.png",
       gg_pcoa_wu_ellipse_arrow,
       height=4, width=6)
