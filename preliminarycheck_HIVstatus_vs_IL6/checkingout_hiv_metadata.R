# Load libraries
library(tidyverse)

# Load file
meta <- "preliminarycheck_HIVstatus_vs_IL6/hiv_metadata.tsv"
hiv_meta <- read_delim(meta, delim="\t")

# Retrieve col names
colnames(hiv_meta)

# Keep only variables needed
hiv_filt <- hiv_meta %>%
  rename(IL_6_pg_mL = `IL-6_pg_mL`) %>% # To avoid problems with ggplot later
  select("sample-id", "original-sample-id", "PID", "HIV_ID",
         "Record_ID", "Cohort", "SEQ", "HIV_Status", 
         "Visit", "Visit_date", "Visit_Age", "Antibiotics", 
         "Blood_Date", "HIV-1_viral_load", "ART_Start_Date", "Current_ART", "IL_6_pg_mL") %>%
  filter(Visit == "2") %>% # Only visit 2 kept due to follow-up inconsistency with visit 3 and stool samples only being collected at visit 2
  arrange(HIV_Status) %>%
  mutate(log_IL6 = log10(IL_6_pg_mL + 1)) # For normality testing

# How many samples are HIV- and HIV+?
for (v in unique(hiv_filt$HIV_Status)) {
  print(v)
  print(nrow(subset(hiv_filt, hiv_filt$HIV_Status == v)))
}

# Plot IL-6 levels of HIV- vs. HIV+
hiv_filt_grouped <- group_by(hiv_filt, HIV_Status)
hivstatus_vs_IL_6 <- ggplot(hiv_filt, aes(x=HIV_Status, y=IL_6_pg_mL)) +
  geom_boxplot() +
  labs(x="HIV status", y="IL-6 levels (pg/mL")

# Normality testing
shapiro.test(hiv_filt$log_IL6)

# Wilcoxon Mann Whitney test
wilcox.test(IL_6_pg_mL ~ HIV_Status, data = hiv_filt)
