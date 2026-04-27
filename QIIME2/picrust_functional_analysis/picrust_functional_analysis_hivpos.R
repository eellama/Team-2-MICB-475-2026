#### Load libraries ####
# install.packages("mlr3tuning")
# install.packages("pheatmap")
library(tidyverse)
library(phyloseq)
library(ggpicrust2)
library("mlr3tuning")
library(ggprism)
library(pheatmap)

load("filtered_and_rarefied/hiv_final.RData")
meta_hivpos <- data.frame(sample_data(hiv_final)) %>%
  rownames_to_column('sample_name') %>%
  mutate(IL6_bin = ifelse(IL.6_pg_mL > 5, "high", "low")) %>%
  filter(HIV_Status == "Positive")

ko <- read.delim("picrust_functional_analysis/pred_metagenome_unstrat.tsv", sep = "\t", row.names = 1)

# Filter out samples in ko that are absent in meta_hivpos
shared_samples_hivpos <- intersect(colnames(ko), meta_hivpos$sample_name)
ko_filt_hivpos <- ko[, shared_samples_hivpos]

#### Differential analysis ####
# Perform pathway differential abundance analysis (DAA) using LinDA method
# Group by IL-6 bin
daa_results_byil6bin_hivpos_df = pathway_daa(abundance = ko,
                                      metadata = meta_hivpos, 
                                      group = "IL6_bin", 
                                      daa_method = "LinDA", 
                                      select = NULL, reference = NULL)
daa_annotated_results_byil6bin_hivpos_df = pathway_annotation(pathway = "KO",
                                                       daa_results_df = daa_results_byil6bin_hivpos_df,
                                                       ko_to_kegg = TRUE)
saveRDS(daa_annotated_results_byil6bin_hivpos_df,'picrust_functional_analysis/daa_annotated_results_byil6bin_hivpos_df.rds')

# Generate pathway error bar plot
daa_annotated_results_byil6bin_hivpos_df_filt <- daa_annotated_results_byil6bin_hivpos_df %>% filter(p_adjust< 0.05) 
nrow(daa_annotated_results_byil6bin_hivpos_df_filt)

peb_il6bin_hivpos = pathway_errorbar_fixed(abundance = ko_filt_hivpos, 
                                     daa_results_df = daa_annotated_results_byil6bin_hivpos_df, 
                                     Group = meta_hivpos$IL6_bin, 
                                     wrap_label = T, wraplength=60,
                                     fc_cutoff = 0, order_by_log = F,
                                     p_values_threshold = 0.05, 
                                     order = "pathway_class", 
                                     ko_to_kegg = FALSE, # changed from TRUE to FALSE so plot can be visualized
                                     p_value_bar = TRUE, 
                                     x_lab = "pathway_name")
peb_il6bin_hivpos

# Save plot
ggsave("picrust_functional_analysis/KEGG_Error_Bar_Fixed_il6bin_hivpos.png",
       plot = peb_il6bin_hivpos, width = 20, height = 10)

# Check for identical features
sig_features_hivpos = daa_annotated_results_byil6bin_hivpos_df %>% 
  pull('feature') %>% unique()

ko_relab = read.delim('picrust_functional_analysis/pred_metagenome_unstrat.tsv',row.names = 1) %>% 
  apply(2,function(x) x/sum(x)) %>% as.data.frame()
colSums(ko_relab) # Should be all 1

stats_spearman_hivpos = ko_relab %>% t() %>%  # switch rows and columns
  as.data.frame() %>% # Lets us use select()
  select(all_of(sig_features_hivpos)) %>% 
  cor(method = 'spearman') 

stats_pearson_hivpos = ko_relab %>% t() %>%  # switch rows and columns
  as.data.frame() %>% # Lets us use select()
  select(all_of(sig_features_hivpos)) %>% 
  cor(method = 'pearson') #Pearson tests between every feature pair

# Define the color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(40)

# Define breaks for the color scale
breaks <- seq(-1, 1, length.out = 41)  # 40 colors + 1 for the endpoint

# Create the clustered heatmap with centered color at zero (can't run this as only 1 sample present)
pheatmap(stats_spearman_hivpos, 
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",          
         color = color_palette, 
         breaks = breaks,
         main = '',
         fontsize_row = 10, 
         fontsize_col = 10,
         filename = "picrust_functional_analysis/KEGG_Heatmap_il6bin_hivpos_spearman.png",
         height= 9, width = 9)

pheatmap(stats_pearson_hivpos, 
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",          
         color = color_palette, 
         breaks = breaks,
         main = '',
         fontsize_row = 10, 
         fontsize_col = 10,
         filename = "picrust_functional_analysis/KEGG_Heatmap_il6bin_hivpos_pearson.png",
         height= 9, width = 9)

# See how similar sample groups are
# Generate pathway PCA plot
Pathway_PCA_il6bin_hivpos <- pathway_pca(abundance = ko_filt_hivpos,
                                   metadata = meta_hivpos, 
                                   group = "IL6_bin")

# Save plot
ggsave("picrust_functional_analysis/Pathway_PCA_il6bin_hivpos.png",
       plot = Pathway_PCA_il6bin_hivpos, width = 10, height = 10)
