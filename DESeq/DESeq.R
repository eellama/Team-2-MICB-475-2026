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
                          contrast = c("HIV_Status","Positive","Negative"))

# Comparing IL6 high vs low
il6_bin_deseq <- phyloseq_to_deseq2(hiv_plus1, ~`IL6_bin`)
DESEQ_il6_bin <- DESeq(il6_bin_deseq)
il6_bin_res <- results(DESEQ_il6_bin, tidy=TRUE, 
                       contrast = c("IL6_bin","high","low"))

# IL6 Hi vs Lo on HIV+ patients
hiv_plwh <- subset_samples(hiv_final, HIV_Status == "Positive")
hiv_plwh_plus1 <- transform_sample_counts(hiv_plwh, function(x) x+1)
hiv_plwh_deseq <- phyloseq_to_deseq2(hiv_plwh_plus1, ~`IL6_bin`)
DESEQ_hiv_plwh <- DESeq(hiv_plwh_deseq)
HIV_plwh_res <- results(DESEQ_hiv_plwh, tidy=TRUE,
                          contrast = c("IL6_bin","high","low"))
# IL6 on healthy patients
hiv_neg <- subset_samples(hiv_final, HIV_Status == "Negative")
hiv_neg_plus1 <- transform_sample_counts(hiv_neg, function(x) x+1)
hiv_neg_deseq <- phyloseq_to_deseq2(hiv_neg_plus1, ~`IL6_bin`)
DESEQ_hiv_neg <- DESeq(hiv_neg_deseq)
HIV_neg_res <- results(DESEQ_hiv_neg, tidy=TRUE,
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

# HIV+ IL6 Binned
HIV_plwh_vol_plot <- HIV_plwh_res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="HIV_plwh_vol_plot.png",HIV_plwh_vol_plot)

# Healthy IL6 Binned
HIV_neg_vol_plot <- HIV_neg_res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="HIV_neg_vol_plot.png",HIV_neg_vol_plot)

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

# PLWH IL-6
HIV_plwh_sigASVs <- HIV_plwh_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(HIV_plwh_sigASVs)
# Get only asv names
HIV_plwh_sigASVs_vec <- HIV_plwh_sigASVs %>%
  pull(ASV)

# Prune phyloseq file
HIV_plwh_deseq <- prune_taxa(HIV_plwh_sigASVs_vec,hiv_final)
HIV_plwh_sigASVs <- tax_table(HIV_plwh_deseq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(HIV_plwh_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

PLWH_bar <- ggplot(HIV_plwh_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="PLWH_bar.png", PLWH_bar)

# Healthy IL-6
HIV_neg_sigASVs <- HIV_neg_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(HIV_neg_sigASVs)
# Get only asv names
HIV_neg_sigASVs_vec <- HIV_neg_sigASVs %>%
  pull(ASV)

# Prune phyloseq file
HIV_neg_deseq <- prune_taxa(HIV_neg_sigASVs_vec,hiv_final)
HIV_neg_sigASVs <- tax_table(HIV_neg_deseq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(HIV_neg_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

healthy_bar <- ggplot(HIV_neg_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="healthy_bar.png", healthy_bar)
