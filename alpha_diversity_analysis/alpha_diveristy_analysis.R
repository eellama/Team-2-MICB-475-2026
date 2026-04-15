#loading required packages
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

#### Load in RData ####
load("hiv_rare.RData")
load("hiv_final.RData")

#Binning IL-6 levels
IL6_binning <- data.frame(sample_data(hiv_rare))
IL6_binned_hiv_rare <- mutate(IL6_binning, IL6_bin = ifelse(IL.6_pg_mL > 5, "high", "low"))

sample_data(hiv_rare) <- sample_data(IL6_binned_hiv_rare)

#filtering for HIV positive and negative
pos_HIV_IL6_binned_rare <- subset_samples(hiv_rare, HIV_Status == "Positive")
neg_HIV_IL6_binned_rare <- subset_samples(hiv_rare, HIV_Status == "Negative")

save(pos_HIV_IL6_binned_rare, file="pos_HIV_IL6_binned_rare.RData")
save(neg_HIV_IL6_binned_rare, file="neg_HIV_IL6_binned_rare.RData")


#### Alpha Diveristy ####
plot_richness(hiv_rare)
plot_richness(hiv_rare, x = "IL6_bin", measures = c("Shannon")) 

gg_richness <- plot_richness(hiv_rare, x = "IL6_bin", measures = c("Shannon")) +
  xlab("IL6 Levels") +
  geom_boxplot() +
  geom_point() +
  ggtitle("Shannon Diversity")
gg_richness

ggsave(filename = "Shannon_Diveristy.png"
       , gg_richness
       , height=4, width=6)


plot_richness(pos_HIV_IL6_binned_rare, x = "IL6_bin", measures = c("Shannon")) +
  xlab("IL6 Levels") +
  geom_boxplot() +
  geom_point() +
  ggtitle("HIV POSITIVE Shannon Diversity")

plot_richness(neg_HIV_IL6_binned_rare, x = "IL6_bin", measures = c("Shannon")) +
  xlab("IL6 Levels") +
  geom_boxplot() +
  geom_point() +
  geom_jitter() +
  ggtitle("HIV NEGATIVE Shannon Diversity")

# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(hiv_rare)), phy_tree(hiv_rare),
                 include.root=F) 
?pd

phylo_dist_pos <- pd(t(otu_table(pos_HIV_IL6_binned_rare)), phy_tree(pos_HIV_IL6_binned_rare),
                       include.root=F) 

phylo_dist_neg <- pd(t(otu_table(neg_HIV_IL6_binned_rare)), phy_tree(neg_HIV_IL6_binned_rare),
                     include.root=F) 

# add PD to metadata table
sample_data(hiv_rare)$PD <- phylo_dist$PD

sample_data(pos_HIV_IL6_binned_rare)$PD <- phylo_dist_pos$PD

sample_data(neg_HIV_IL6_binned_rare)$PD <- phylo_dist_neg$PD

# calculate Shannon's Diversity 
estimate_richness(hiv_rare, measures = "Shannon")
estimate_richness(pos_HIV_IL6_binned_rare, measures = "Shannon")
estimate_richness(neg_HIV_IL6_binned_rare, measures = "Shannon")

shannon_d <- estimate_richness(hiv_rare, measures = "Shannon")
shannon_d_pos <- estimate_richness(pos_HIV_IL6_binned_rare, measures = "Shannon")
shannon_d_neg <- estimate_richness(neg_HIV_IL6_binned_rare, measures = "Shannon")

sample_data(hiv_rare)$Shannon <- shannon_d$Shannon
sample_data(pos_HIV_IL6_binned_rare)$Shannon <- shannon_d_pos$Shannon
sample_data(neg_HIV_IL6_binned_rare)$Shannon <- shannon_d_neg$Shannon




# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(hiv_rare), aes(IL6_bin, PD)) + 
  geom_boxplot() +
  geom_point() +
  xlab("IL6 Levels") +
  ylab("Phylogenetic Diversity") +
  ggtitle("Faith's PD")

# view plot
plot.pd

#saving plot
ggsave(filename = "FaithsPD.png", plot.pd)

#positive HIV
plot.pdpos <- ggplot(sample_data(pos_HIV_IL6_binned_rare), aes(IL6_bin, PD)) + 
  geom_boxplot() +
  geom_point() +
  xlab("IL6 Levels") +
  ylab("Phylogenetic Diversity") +
  ggtitle("HIV POSITIVE Faith's PD")
plot.pdpos

#negative HIV
plot.pdneg <- ggplot(sample_data(neg_HIV_IL6_binned_rare), aes(IL6_bin, PD)) + 
  geom_boxplot() +
  geom_point() +
  xlab("IL6 Levels") +
  ylab("Phylogenetic Diversity") +
  ggtitle("HIV NEGATIVE Faith's PD")
plot.pdneg

#### Stats ####
samp_dat <- data.frame(sample_data(hiv_rare))
samp_datpos <- data.frame(sample_data(pos_HIV_IL6_binned_rare))
samp_datneg <- data.frame(sample_data(neg_HIV_IL6_binned_rare))

#Kruskal-Wallis rank sum test
kruskal.test(PD ~ IL6_bin, data=samp_dat)

#positive HIV
kruskal.test(PD ~ IL6_bin, data = samp_datpos)

#negative HIV
kruskal.test(PD ~ IL6_bin, data = samp_datneg)


#Spearman rank correlation
cor.test(samp_dat$IL.6_pg_mL, samp_dat$PD, method="spearman")
cor.test(samp_datpos$IL.6_pg_mL, samp_datpos$PD, method="spearman")
cor.test(samp_datneg$IL.6_pg_mL, samp_datneg$PD, method="spearman")
cor.test(samp_datpos$IL.6_pg_mL, samp_datpos$Shannon, method="spearman")
cor.test(samp_datneg$IL.6_pg_mL, samp_datneg$Shannon, method="spearman")

#PD vs IL6 levels plots
#negative 
plot.pdnegss <- ggplot(sample_data(neg_HIV_IL6_binned_rare), aes(IL.6_pg_mL, PD)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("IL6 Levels") +
  ylab("Phylogenetic Diversity") +
  theme_classic() +
  ggtitle("HIV NEGATIVE - Faith's PD")
plot.pdnegss
ggsave(filename = "HIV_negative - FaithsPDvsIL6.png"
       , plot.pdnegss
       , height=4, width=6)
#positive
plot.pdposss <- ggplot(sample_data(pos_HIV_IL6_binned_rare), aes(IL.6_pg_mL, PD)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("IL6 Levels") +
  ylab("Phylogenetic Diversity") +
  theme_classic() +
  ggtitle("HIV POSITIVE - Faith's PD")
plot.pdposss
ggsave(filename = "HIV_positive - FaithsPDvsIL6.png"
       , plot.pdposss
       , height=4, width=6)
#

#all
plot.pdall <- ggplot(sample_data(hiv_rare), aes(IL.6_pg_mL, PD)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab ("IL6 Levels") +
  ylab("Phylogenetic Diversity") +
  theme_classic() +
  ggtitle ("All samples - Faith's PD")
plot.pdall
ggsave(filename = "All_samples - FaithsPDvsIL6.png"
       , plot.pdall
       , height=4, width=6)
#
#Shannon vs IL6 levels plots
#negative
plot_shannon_neg <- ggplot(sample_data(neg_HIV_IL6_binned_rare), aes(IL.6_pg_mL, Shannon)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("IL6 Levels") +
  ylab("Shannon's Diversity") +
  theme_classic() +
  ggtitle("HIV NEGATIVE - Shannon's Diversity")
plot_shannon_neg
ggsave(filename = "HIV_negative - ShannonvsIL6.png"
       , plot_shannon_neg
       , height=4, width=6)
#negative
plot_shannon_pos <- ggplot(sample_data(pos_HIV_IL6_binned_rare), aes(IL.6_pg_mL, Shannon)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("IL6 Levels") +
  ylab("Shannon's Diversity") +
  theme_classic() +
  ggtitle("HIV POSITIVE - Shannon's Diversity")
plot_shannon_pos
ggsave(filename = "HIV_positive - ShannonsvsIL6.png"
       , plot_shannon_pos
       , height=4, width=6)
#all
plot_shannon_all <- ggplot(sample_data(hiv_rare), aes(IL.6_pg_mL, Shannon)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("IL6 Levels") +
  ylab("Shannon's Diversity") +
  theme_classic() +
  ggtitle("All samples - Shannon's Diversity")
plot_shannon_all
ggsave(filename = "All_samples - ShanonsvsIL6.png"
       , plot_shannon_all
       , height=4, width=6)

#Pearson correlation
cor.test(samp_dat$IL.6_pg_mL, samp_dat$PD, method="pearson")
cor.test(samp_datpos$IL.6_pg_mL, samp_datpos$PD, method="pearson")
cor.test(samp_datneg$IL.6_pg_mL, samp_datneg$PD, method="pearson")

#other stats for fun
#pearson correlation
cor.test(samp_dat$IL.6_pg_mL, samp_dat$PD, method="pearson")
#Wilcoxon Rank Sum Test for fun
wilcox.test(PD ~ IL6_bin, data=samp_dat, exact = FALSE)
#plotting IL6 against alpha diversity (faith's)
ggplot(samp_dat,aes(x=IL.6_pg_mL, y=PD)) +
  geom_point() +
  geom_smooth(method = "lm")

