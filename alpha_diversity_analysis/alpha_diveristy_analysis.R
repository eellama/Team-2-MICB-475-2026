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

# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(hiv_rare)), phy_tree(hiv_rare),
                 include.root=F) 
?pd
# add PD to metadata table
sample_data(hiv_rare)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(hiv_rare), aes(IL6_bin, PD)) + 
  geom_boxplot() +
  geom_point() +
  xlab("IL6 Levels") +
  ylab("Phylogenetic Diversity") +
  ggtitle("Faith's PD")

# view plot
plot.pd

ggsave(filename = "FaithsPD.png", plot.pd)

#### Stats ####
samp_dat <- data.frame(sample_data(hiv_rare))
#Kruskal-Wallis rank sum test
kruskal.test(PD ~ IL6_bin, data=samp_dat)


#Spearman rank correlation
cor.test(samp_dat$IL.6_pg_mL, samp_dat$PD, method="spearman")




#other stats for fun
#pearson correlation
cor.test(samp_dat$IL.6_pg_mL, samp_dat$PD, method="pearson")
#Wilcoxon Rank Sum Test for fun
wilcox.test(PD ~ IL6_bin, data=samp_dat, exact = FALSE)
#plotting IL6 against alpha diversity (faith's)
ggplot(samp_dat,aes(x=IL.6_pg_mL, y=PD)) +
  geom_point() +
  geom_smooth(method = "lm")

