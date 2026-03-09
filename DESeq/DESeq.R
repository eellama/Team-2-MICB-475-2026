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

# Just comparing HIV_Status
hiv_plus1 <- transform_sample_counts(hiv_final, function(x) x+1)
hiv_status_deseq <- phyloseq_to_deseq2(hiv_plus1, ~`HIV_Status`)
DESEQ_hiv_status <- DESeq(hiv_status_deseq)
HIV_status_res <- results(DESEQ_hiv_status, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("HIV_Status","Positive","Negative"))

# comparing IL6 high vs low
il6_bin_deseq <- phyloseq_to_deseq2(hiv_plus1, ~`IL6_bin`)
DESEQ_il6_bin <- DESeq(il6_bin_deseq)
il6_bin_res <- results(DESEQ_il6_bin, tidy=TRUE, 
                          #this will ensure that No is your reference group
                          contrast = c("IL6_bin","high","low"))

#### Visualizing Data ####

## Volcano plot: effect size VS significance

# HIV_Status
HIV_status_vol_plot <- HIV_status_res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="HIV_status_vol_plot.png",HIV_status_vol_plot)

# IL6 Binned
il6_bin_vol_plot <- il6_bin_res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="il6_bin_vol_plot.png",il6_bin_vol_plot)

## Bar graph

# HIV_Status
HIV_status_sigASVs <- HIV_status_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(HIV_status_sigASVs)
# Get only asv names
HIV_status_sigASVs_vec <- HIV_status_sigASVs %>%
  pull(ASV)

# Prune phyloseq file
hiv_status_deseq <- prune_taxa(HIV_status_sigASVs_vec,hiv_final)
HIV_status_sigASVs <- tax_table(hiv_status_deseq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(HIV_status_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

HIV_status_bar <- ggplot(HIV_status_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="HIV_status_bar.png", HIV_status_bar)

# IL6 Binned
il6_bin_sigASVs <- il6_bin_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(il6_bin_sigASVs)
# Get only asv names
il6_bin_sigASVs_vec <- il6_bin_sigASVs %>%
  pull(ASV)

# Prune phyloseq file
il6_bin_deseq <- prune_taxa(il6_bin_sigASVs_vec,hiv_final)
il6_bin_sigASVs <- tax_table(il6_bin_deseq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(il6_bin_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

IL6_bar <- ggplot(il6_bin_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="IL6_bar.png", IL6_bar)
