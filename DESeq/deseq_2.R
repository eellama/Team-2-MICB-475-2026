#### Load packages ####
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(ggplot2)

#### Load data ####
load("hiv_final.RData")

IL6_binning_update <- data.frame(sample_data(hiv_final)) %>%
  mutate(IL6_bin = ifelse(IL.6_pg_mL > 5, "high", "low"))
sample_data(hiv_final) <- sample_data(IL6_binning_update)

#### DESeq ####

# IL6 Hi vs Lo on HIV+ patients
hiv_plwh <- subset_samples(hiv_final, HIV_Status == "Positive")
hiv_plwh_plus1 <- transform_sample_counts(hiv_plwh, function(x) x+1)
hiv_plwh_deseq <- phyloseq_to_deseq2(hiv_plwh_plus1, ~`IL6_bin`)
DESEQ_hiv_plwh <- DESeq(hiv_plwh_deseq)
HIV_plwh_res <- results(DESEQ_hiv_plwh, tidy=TRUE,
                        contrast = c("IL6_bin","high","low"))

# glom to Genus
HIV_plwh_ASVs <- HIV_plwh_res %>%
  dplyr::rename(ASV=row)
HIV_plwh_res_vec <- HIV_plwh_ASVs %>%
  pull(ASV)
HIV_plwh_genus_deseq <- prune_taxa(HIV_plwh_res_vec,hiv_final)
HIV_plwh_ASVs <- tax_table(HIV_plwh_genus_deseq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(HIV_plwh_ASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

#### Visualizing Data ####
# Prune phyloseq file
HIV_plwh_sigASVs <- HIV_plwh_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(HIV_plwh_sigASVs)
# Get only significant asv names
HIV_plwh_sigASVs_vec <- HIV_plwh_sigASVs %>%
  pull(ASV)

HIV_plwh_deseq <- prune_taxa(HIV_plwh_sigASVs_vec,hiv_final)
HIV_plwh_sigASVs <- tax_table(HIV_plwh_deseq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(HIV_plwh_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

## Volcano plot: effect size VS significance
# HIV+ IL6 Binned
HIV_plwh_genus_vol_plot <- HIV_plwh_ASVs %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="HIV_plwh_genus_vol_plot.png",HIV_plwh_genus_vol_plot)

## Bar graph
# PLWH IL-6
PLWH_bar <- ggplot(HIV_plwh_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="PLWH_bar.png", PLWH_bar, width = 30, height = 15, units = "cm")
