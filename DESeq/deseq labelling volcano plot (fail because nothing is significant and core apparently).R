#### Load packages ####
library(patchwork)
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(microbiome)
library(ggrepel)

#### Load data ####
load("hiv_final.RData")

IL6_binning_update <- data.frame(sample_data(hiv_final)) %>%
  mutate(IL6_bin = ifelse(IL.6_pg_mL > 5, "high", "low"))
sample_data(hiv_final) <- sample_data(IL6_binning_update)

#### Core ####
# subset IL-6 bins for core microbiome
il6_low <- subset_samples(hiv_final, `IL6_bin`=="low")
il6_high <- subset_samples(hiv_final, `IL6_bin`=="high")

# identify members in each group that meet thresholds
low_ASVs <- core_members(il6_low, detection=0, prevalence = 0.7)
high_ASVs <- core_members(il6_high, detection=0, prevalence = 0.7)

ggVennDiagram(x=list(high_ASVs, low_ASVs))

# Get taxonomy from full dataset
tax_df_full <- tax_table(hiv_final) %>%
  as.data.frame() %>%
  rownames_to_column("ASV")

# Map ASVs to Genus
low_core_genus <- tax_df_full %>%
  filter(ASV %in% low_ASVs) %>%
  pull(Genus) %>%
  unique()

high_core_genus <- tax_df_full %>%
  filter(ASV %in% high_ASVs) %>%
  pull(Genus) %>%
  unique()

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

HIV_plwh_res <- HIV_plwh_res %>%
  mutate(core_status = case_when(
    Genus %in% high_core_genus & Genus %in% low_core_genus ~ "Core Both",
    Genus %in% high_core_genus ~ "Core High",
    Genus %in% low_core_genus ~ "Core Low",
    TRUE ~ "Non-core"
  ))

#### Visualizing Data ####

## Volcano plot: effect size VS significance
# HIV+ IL6 Binned
HIV_plwh_genus_vol <- HIV_plwh_res %>%
  mutate(category = case_when(log2FoldChange > 2  ~ "Upregulated",
                              log2FoldChange < -2 ~ "Downregulated",
                              padj < 0.05 ~ "Significant",
                              TRUE ~ "Not Significant"))

genes_to_label <- HIV_plwh_genus_vol[HIV_plwh_genus_vol$core_status != "Non-core", "Genus"]

HIV_plwh_genus_vol_plot <- ggplot(HIV_plwh_genus_vol,
                                  aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = core_status)) +
  ggrepel::geom_label_repel(data = subset(HIV_plwh_genus_vol, padj < 0.05 & abs(log2FoldChange) > 2 & core_status != "Non-core"), 
                            aes(label = Genus), 
                            color = "black", 
                            show.legend = FALSE, 
                            max.overlaps = Inf) +
  theme_bw()

ggsave(filename="HIV_plwh_genus_vol_plot.png",HIV_plwh_genus_vol_plot)

## Bar graph
# Filter for significant genera
HIV_plwh_sig_genus <- HIV_plwh_res %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2)

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

HIV_res <- HIV_res %>%
  mutate(core_status = case_when(
    Genus %in% high_core_genus & Genus %in% low_core_genus ~ "Core Both",
    Genus %in% high_core_genus ~ "Core High",
    Genus %in% low_core_genus ~ "Core Low",
    TRUE ~ "Non-core"
  ))

#### Visualizing Data ####

## Volcano plot: effect size VS significance
# HIV+ IL6 Binned

HIV_vol <- HIV_res %>%
  mutate(category = case_when(log2FoldChange > 2  ~ "Upregulated",
                              log2FoldChange < -2 ~ "Downregulated",
                              padj < 0.05 ~ "Significant",
                              TRUE ~ "Not Significant"))

genes_label <- HIV_vol[HIV_vol$core_status != "Non-core", "Genus"]

HIV_vol_plot <- ggplot(HIV_vol,
                       aes(x = log2FoldChange, y = -log10(padj))) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey") +
  
  geom_point(aes(color = core_status)) + # in theory I would have it coloured by significance and just label the core but there is no core
  ggrepel::geom_label_repel(data = subset(HIV_vol, padj < 0.05 & abs(log2FoldChange) > 2 & core_status != "Non-core"), 
                            aes(label = Genus), 
                            color = "black", 
                            show.legend = FALSE, 
                            max.overlaps = Inf) +
  theme_bw()

ggsave(filename="all_status_vol_plot.png",HIV_vol_plot)

## Bar graph
# Filter for significant genera
HIV_sig_genus <- HIV_res %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2)

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
volcanos <- HIV_vol_plot + HIV_plwh_genus_vol_plot

ggsave(filename="bars.png", bars)
ggsave(filename="volcanos.png", volcanos)

