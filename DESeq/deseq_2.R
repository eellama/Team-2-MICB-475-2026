#### Load packages ####

library(patchwork)
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
hiv_plwh_genus <- tax_glom(hiv_plwh_plus1, taxrank = "Genus")
hiv_plwh_deseq <- phyloseq_to_deseq2(hiv_plwh_genus, ~`IL6_bin`)
DESEQ_hiv_plwh <- DESeq(hiv_plwh_deseq)
HIV_plwh_res <- results(DESEQ_hiv_plwh, tidy=TRUE,
                        contrast = c("IL6_bin","high","low"))

# Glom to Genus
HIV_plwh_res <- HIV_plwh_res %>%
  dplyr::rename(Genus_ID=row)

tax_df <- tax_table(hiv_plwh_genus) %>%
  as.data.frame() %>%
  rownames_to_column("Genus_ID")

HIV_plwh_res <- right_join(HIV_plwh_res, tax_df, by = "Genus_ID") %>%
  mutate(Genus = make.unique(as.character(Genus)))

#### Visualizing Data ####

## Volcano plot: effect size VS significance
# HIV+ IL6 Binned
HIV_plwh_genus_vol_plot <- HIV_plwh_res %>%
  mutate(category = case_when(padj < 0.05 & log2FoldChange > 1.5  ~ "Upregulated",
                              padj < 0.05 & log2FoldChange < -1.5 ~ "Downregulated",
                              TRUE ~ "Not Significant")) %>%
  ggplot() +
  geom_vline(xintercept = -1.5, linetype = "dashed", colour = "grey") +
  geom_vline(xintercept = 1.5, linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey") +
  scale_x_continuous(breaks = c(-3, -1.5, 0, 1.5)) +
  labs(x = "Fold Change (log[2])", y = "-log[10] (Padj)") +
  scale_color_manual(
    values = c("Upregulated" = "green",
               "Downregulated" = "#e31a1c",
               "Not Significant" = "#1f78b4"), 
    name = "") +
  geom_point(aes(x=log2FoldChange, y= -log10(padj), col=category)) +
  theme_bw(9) # remove background gray grid
ggsave(filename="HIV_plwh_genus_vol_plot.png",HIV_plwh_genus_vol_plot)

## Bar graph
# Filter for significant genera
HIV_plwh_sig_genus <- HIV_plwh_res %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5)

# PLWH IL-6
PLWH_bar <- ggplot(HIV_plwh_sig_genus) +
  geom_bar(aes(x=reorder(Genus, log2FoldChange), y=log2FoldChange, fill = log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  theme_bw() +
  labs(x = "Genus", y = "log[2] Fold Change")
ggsave(filename="PLWH_bar.png", PLWH_bar)

#### ALL STATUS ####
#### DESeq ####

# IL6 Hi vs Lo on all patients
hiv_plus1 <- transform_sample_counts(hiv_final, function(x) x+1)
hiv_genus <- tax_glom(hiv_plus1, taxrank = "Genus")
hiv_deseq <- phyloseq_to_deseq2(hiv_genus, ~`IL6_bin`)
DESEQ_hiv <- DESeq(hiv_deseq)
HIV_res <- results(DESEQ_hiv, tidy=TRUE,
                        contrast = c("IL6_bin","high","low"))

# Glom to Genus
HIV_res <- HIV_res %>%
  dplyr::rename(Genus_ID=row)

tax_all_df <- tax_table(hiv_genus) %>%
  as.data.frame() %>%
  rownames_to_column("Genus_ID")

HIV_res <- right_join(HIV_res, tax_all_df, by = "Genus_ID") %>%
  mutate(Genus = make.unique(as.character(Genus)))

#### Visualizing Data ####

## Volcano plot: effect size VS significance
# HIV+ IL6 Binned
HIV_vol_plot <- HIV_res %>%
  mutate(category = case_when(padj < 0.05 & log2FoldChange > 1.5  ~ "Upregulated",
                              padj < 0.05 & log2FoldChange < -1.5 ~ "Downregulated",
                              TRUE ~ "Not Significant")) %>%
  ggplot() +
  geom_vline(xintercept = -1.5, linetype = "dashed", colour = "grey") +
  geom_vline(xintercept = 1.5, linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey") +
  scale_x_continuous(breaks = c(-3, -1.5, 0, 1.5)) +
  labs(x = "Fold Change (log[2])", y = "-log[10] (Padj)") +
  geom_point(aes(x=log2FoldChange, y= -log10(padj), col=category)) +
  scale_color_manual(
    values = c("Upregulated" = "green",
               "Downregulated" = "#e31a1c",
               "Not Significant" = "#1f78b4"), 
    name = "") +
  theme_bw(9) # remove background gray grid
ggsave(filename="all_status_vol_plot.png",HIV_vol_plot)

## Bar graph
# Filter for significant genera
HIV_sig_genus <- HIV_res %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5)

# PLWH IL-6
all_bar <- ggplot(HIV_sig_genus) +
  geom_bar(aes(x=reorder(Genus, log2FoldChange), y=log2FoldChange, fill = log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  theme_bw() +
  labs(x = "Genus", y = "log[2] Fold Change")
ggsave(filename="all_bar.png", all_bar)

# Put together
bars <- all_bar + PLWH_bar
volcanos <- HIV_vol_plot / HIV_plwh_genus_vol_plot

ggsave(filename="bars.png", bars, width = 3000, height = 2500, units = "px")
ggsave(filename="volcanos.png", volcanos, width = 3000, height = 2500, units = "px")
