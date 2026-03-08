library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

load("hiv_final.RData")
# This file was uploaded on Feb25 2026 onto the group github. 

### Core Microbiome #1: Simple IL6 high vs IL6 low, all persons included.

#convert the file to relative abundance. 
hiv_RA <- transform_sample_counts(hiv_final, fun=function(x) x/sum(x))

#bin person by IL-6 level: <5pg/mL or >5pg/mL
#included "sample_data" command here because I was having issues with making this command
#actually bin stuff. 
lowil6 <- subset_samples(hiv_RA, sample_data(hiv_RA)$'IL.6_pg_mL' <= 5)
highil6 <- subset_samples(hiv_RA, sample_data(hiv_RA)$'IL.6_pg_mL' >= 5)

#Which ASVs are in >70% of samples of the subsets?
low_ASVs <- core_members(lowil6, detection=0, prevalence = 0.7)
high_ASVs <- core_members(highil6, detection=0, prevalence = 0.7)

#Venn Diagram
il_6list <- list(IL6high = high_ASVs, IL6low = low_ASVs)
allpatientvenn1 <- ggVennDiagram(x = il_6list) + ggtitle("Core ASVs in IL6 High vs Low groups\n(70% prevalence), all patients")
ggsave("allpatientvenn1.jpeg", allpatientvenn1)

#Which ASVs are in >50% of samples of the subsets? Wanted to see just out of curiousity
low_ASVs <- core_members(lowil6, detection=0, prevalence = 0.5)
high_ASVs <- core_members(highil6, detection=0, prevalence = 0.5)

il_6list <- list(IL6high = high_ASVs, IL6low = low_ASVs)
allpatientvenn2 <- ggVennDiagram(x = il_6list) + ggtitle("Core ASVs in IL6 High vs Low groups\n(50% prevalence), all patients")
ggsave("allpatientvenn2.jpeg", allpatientvenn2)




### Core Microbiome #2:IL-6 high v IL-6 low in HIV+ persons ONLY
hiv_only <- subset_samples(mpt_RA, HIV_diagnosis == "Negative") #only keeps HIV+ 
HPlowil6 <- subset_samples(hiv_only, sample_data(hiv_only)$'IL.6_pg_mL' <= 5) #perform same IL6 binning as previous. 
HPhighil6 <- subset_samples(hiv_only, sample_data(hiv_only)$'IL.6_pg_mL' >= 5)

#Which ASVs are in >70% of samples of the subsets?
HPlow_ASVs <- core_members(HPlowil6, detection=0, prevalence = 0.7)
HPhigh_ASVs <- core_members(HPhighil6, detection=0, prevalence = 0.7)

#Venn Diagram
il_6list <- list(IL6high = HPhigh_ASVs, IL6low = HPlow_ASVs)
HPvenn1 <- ggVennDiagram(x = il_6list) + ggtitle("Core ASVs in IL6 High vs Low groups\n(70% prevalence), HIV+ only")
ggsave("HPvenn1.jpeg", HPvenn1)

#Which ASVs are in >50% of samples of the subsets? Wanted to see just out of curiousity
HPlow_ASVs <- core_members(HPlowil6, detection=0, prevalence = 0.5)
HPhigh_ASVs <- core_members(HPhighil6, detection=0, prevalence = 0.5)

il_6list <- list(IL6high = HPhigh_ASVs, IL6low = HPlow_ASVs)
HPvenn2 <- ggVennDiagram(x = il_6list) + ggtitle("Core ASVs in IL6 High vs Low groups\n(50% prevalence), HIV+ only")
ggsave("HPvenn2.jpeg", HPvenn2)

