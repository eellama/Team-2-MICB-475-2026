#### Load libraries ####
# install.packages("mlr3tuning")
# install.packages("pheatmap")
library(tidyverse)
library(phyloseq)
library(ggpicrust2)
library("mlr3tuning")
library(ggprism)
library(pheatmap)

#### Load data ####
load("filtered_and_rarefied/hiv_final.RData")
meta <- data.frame(sample_data(hiv_final)) %>%
  rownames_to_column('sample_name') %>%
  mutate(IL6_bin = ifelse(IL.6_pg_mL > 5, "high", "low"))

ko <- read.delim("picrust_functional_analysis/pred_metagenome_unstrat.tsv", sep = "\t", row.names = 1)

# Filter out samples in ko that are absent in meta
shared_samples <- intersect(colnames(ko), meta$sample_name)
ko_filt <- ko[, shared_samples]

#### Differential analysis ####
# Perform pathway differential abundance analysis (DAA) using LinDA method
# Group by HIV status (annotated results only yield 14 rows. Why?)
daa_results_byhivstat_df = pathway_daa(abundance = ko,
                             metadata = meta, 
                             group = "HIV_Status", 
                             daa_method = "LinDA", 
                             select = NULL, reference = NULL)
daa_annotated_results_byhivstat_df = pathway_annotation(pathway = "KO",
                                              daa_results_df = daa_results_byhivstat_df,
                                              ko_to_kegg = TRUE)
saveRDS(daa_annotated_results_byhivstat_df,'picrust_functional_analysis/daa_annotated_results_byhivstat_df.rds')

# Group by IL-6 bin
daa_results_byil6bin_df = pathway_daa(abundance = ko,
                             metadata = meta, 
                             group = "IL6_bin", 
                             daa_method = "LinDA", 
                             select = NULL, reference = NULL)
daa_annotated_results_byil6bin_df = pathway_annotation(pathway = "KO",
                                              daa_results_df = daa_results_byil6bin_df,
                                              ko_to_kegg = TRUE) # pathway details are NA!
saveRDS(daa_annotated_results_byil6bin_df,'picrust_functional_analysis/daa_annotated_results_byil6bin_df.rds')

# Load .RDS files
daa_annotated_results_byhivstat_df = readRDS('picrust_functional_analysis/daa_annotated_results_byhivstat_df.rds')
daa_annotated_results_byil6bin_df = readRDS('picrust_functional_analysis/daa_annotated_results_byil6bin_df.rds')

# Generate pathway error bar plot
source('picrust_functional_analysis/ggpicrust2_errorbar_function_fixed.R')

daa_annotated_results_byhivstat_df <- daa_annotated_results_byhivstat_df %>% 
  filter(p_adjust< 0.05, 
         abs(log2_fold_change) > 2) # Can choose to add abs(log2_fold_change) > 2 as a second filter
nrow(daa_annotated_results_byhivstat_df)
# daa_annotated_results_byil6bin_df_filt <- daa_annotated_results_byil6bin_df %>% filter(p_adjust< 0.05) 
# nrow(daa_annotated_results_byil6bin_df_filt) # no rows match filter conditions :(

peb_hivstat = pathway_errorbar_fixed(abundance = ko_filt, 
                             daa_results_df = daa_annotated_results_byhivstat_df, 
                             Group = meta$HIV_Status, 
                             wrap_label = T, wraplength=60,
                             fc_cutoff = 0, order_by_log = F,
                             p_values_threshold = 0.05, 
                             order = "pathway_class", 
                             ko_to_kegg = FALSE, # changed from TRUE to FALSE so plot can be visualized
                             p_value_bar = TRUE, 
                             x_lab = "pathway_name")
peb_hivstat

# Save plot
ggsave("picrust_functional_analysis/KEGG_Error_Bar_Fixed_hivstat.png",
       plot = peb_hivstat, width = 20, height = 10)

# Check for identical features
sig_features = daa_annotated_results_byhivstat_df %>% 
  pull('feature') %>% unique()

ko_relab = read.delim('picrust_functional_analysis/pred_metagenome_unstrat.tsv',row.names = 1) %>% 
  apply(2,function(x) x/sum(x)) %>% as.data.frame()
colSums(ko_relab) # Should be all 1

stats_spearman = ko_relab %>% t() %>%  # switch rows and columns
  as.data.frame() %>% # Lets us use select()
  select(all_of(sig_features)) %>% 
  cor(method = 'spearman') 

stats_pearson = ko_relab %>% t() %>%  # switch rows and columns
  as.data.frame() %>% # Lets us use select()
  select(all_of(sig_features)) %>% 
  cor(method = 'pearson') #Pearson tests between every feature pair

# Define the color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(40)

# Define breaks for the color scale
breaks <- seq(-1, 1, length.out = 41)  # 40 colors + 1 for the endpoint

# Create the clustered heatmap with centered color at zero
pheatmap(stats_spearman, 
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",          
         color = color_palette, 
         breaks = breaks,
         main = '',
         fontsize_row = 10, 
         fontsize_col = 10,
         filename = "picrust_functional_analysis/KEGG_Heatmap_hivstat_spearman.png",
         height= 9, width = 9)

pheatmap(stats_pearson, 
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",          
         color = color_palette, 
         breaks = breaks,
         main = '',
         fontsize_row = 10, 
         fontsize_col = 10,
         filename = "picrust_functional_analysis/KEGG_Heatmap_hivstat_pearson.png",
         height= 9, width = 9)

# See how similar your sample groups are. Unfortunately, the size of the text can't be increased, even with theme() commands.
# Generate pathway PCA plot
Pathway_PCA_hivstat <- pathway_pca(abundance = ko_filt,
            metadata = meta, 
            group = "HIV_Status")

# Save plot
ggsave("picrust_functional_analysis/Pathway_PCA_hivstat.png",
       plot = Pathway_PCA_hivstat, width = 10, height = 10)
