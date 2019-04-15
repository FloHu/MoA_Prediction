# aim: get MIC data for all drugs
# based on Chemicals used-Ecoli_MICs.xls

# setup -----------------
rm(list = ls())
library(tidyverse)
library(readxl)
source("./R/misc.R")

# drug_synonyms_file <- "/Volumes/typas/Florian/dbsetup_tables/drug_synonyms.csv"
# drug_synonyms <- read_delim(drug_synonyms_file, delim = ";")
drug_synonyms <- read_delim("/Volumes/typas/Florian/dbsetup_tables_new/drug_synonyms.csv", 
  delim = ";")
mics <- read_excel("./data/Chemicals used-Ecoli_MICs.xls", range = "A1:L112")

# mic_tested means that mics were tested again and if they deviated were recorded
mics <- mics[, c("Antibiotic/Condition", "MIC", "Tested", "...12")]
mics <- rename(mics, 'drugname_typaslab' = 'Antibiotic/Condition', 'mic' = 'MIC', 
  'mic_tested' = 'Tested', 'comments' = '...12')

# first problem: map the drugnames to drugname_typaslab

new_synonyms <- c('Bile acids' = 'BILE', 'Triton X-100' = 'TRITONX', 
  'CHIR- 090' = 'CHIR090', 'Epigallocatechin gallate (EGCG)' = 'EGCG', 
  'EPINEPHRINE' = 'EPINEPHRINE', 'Epinephrine' = 'EPINEPHRINE', 
  'NOREPINEPHRINE' = 'NOREPINEPHRINE', 'Norepinephrine' = 'NOREPINEPHRINE', 
  'EtOH' = 'ETHANOL', 'Nitrofurnatoin' = 'NITROFURANTOIN', 
  'Paraquat dichloride' = 'PARAQUAT', 'Phenazine methosulfate (PMS)' = 'PMS', 
  'CCCP (Carbonyl cyanide 3-chlorophenylhydrazone)' = 'CCCP', 
  'Triclosan/Irgasan' = 'TRICLOSAN', 'Irgasan' = 'TRICLOSAN', 
  'Hydrogen peroxide' = 'PEROXIDE', 'Cefotaxime' = 'CEFOTAXIME', 
  'Fosfomycin' = 'FOSFOMYCIN', 'Erythromycin' = 'ERYTHROMYCIN', 
  'Clarithromycin' = 'CLARITHROMYCIN', 'Nalidixic acid' = 'NALIDIXICACID', 
  'Polymyxin B' = 'POLYMYXINB', 'Chlorpromazine' = 'CHLORPROMAZINE', 
  'Fusidic acid' = 'FUSIDICACID')

new_synonyms <- data.frame(drug_synonym = names(new_synonyms), 
  drugname_typaslab = new_synonyms, stringsAsFactors = FALSE)

drug_synonyms <- rbind(drug_synonyms, new_synonyms)

# exclude cases that don't appear:
anti_join(mics, drug_synonyms, by = c('drugname_typaslab' = 'drug_synonym'))

mics <- semi_join(mics, drug_synonyms, by = c('drugname_typaslab' = 'drug_synonym'))
mics$drugname_typaslab <- 
  drug_synonyms$drugname_typaslab[match(x = mics$drugname_typaslab, 
    table = drug_synonyms$drug_synonym)]

mics$resistant <- ifelse(mics$mic == "resistant", TRUE, FALSE)
mics$mic_parsed <- parse_number(mics$mic, na = c("?", "-", "NA", "resistant"))
mics$mic_tested_parsed <- parse_number(mics$mic_tested, na = c("-", "NA", "yes"))
mics$mic_curated <- mics$mic_parsed
mics$mic_curated[!is.na(mics$mic_tested_parsed)] <- mics$mic_tested_parsed[!is.na(mics$mic_tested_parsed)]
mics <- select(mics, drugname_typaslab, mic, mic_parsed, mic_tested, mic_tested_parsed, 
  mic_curated, resistant, comments)

arrange(mics, drugname_typaslab) %>% print(n = 100)

corrections <- c("A22" = 2, "TRICLOSAN" = 0.8, "EPINEPHRINE" = NA, 
  "NOREPINEPHRINE" = NA, "NOVOBIOCIN" = 40, "CYCLOSERINED" = 10, 
  "POLYMYXINB" = 0.5, "BACITRACIN" = NA, "STREPTOMYCIN" = NA, "SDS" = 5, 
  "BENZALKONIUM" = 20, "CCCP" = 50, "FUSIDICACID" = 200, "STREPTOMYCIN" = 8)

mics$mic_curated[match(names(corrections), mics$drugname_typaslab)] <- corrections
mics <- select(mics, drugname_typaslab, mic_curated, resistant, comments)

write_delim(mics, path = "./data/programmatic_output/MICs.csv", delim = ";")

