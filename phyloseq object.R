# This is just for Elena's computer that is stuck in Chinese and I cannot read the console
# Sys.setenv(LANG="en")

# Load packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#### Load data ####
metafp <- "hiv_metadata.csv"
meta <- read_delim(metafp, delim=",")

otufp <- "export_HIV/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "export_HIV/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "export_HIV/rootedtree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep=";"
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
mpt <- phyloseq(OTU, SAMP, TAX, phylotree)