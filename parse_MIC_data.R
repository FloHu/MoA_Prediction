# aim: get MIC data for all drugs
# based on Chemicals used-Ecoli_MICs.xls

# setup -----------------
rm(list = ls())
library(tidyverse)
library(readxl)
source("./R/myhead.R")

drug_synonyms_file <- "/Volumes/typas/Florian/dbsetup_tables/drug_synonyms.csv"
nichols_mic_file <- "./data/Chemicals used-Ecoli_MICs.xls"
drug_synonyms <- read_delim(drug_synonyms_file, delim = ";")

# mic_tested means that mics were tested again and if they deviated were recorded
mics <- read_excel(nichols_mic_file, range = "A1:L112")
mics <- mics[, c("Antibiotic/Condition", "MIC", "Tested", "X__1")]
mics <- rename(mics, 'drugname_typaslab' = 'Antibiotic/Condition', 'mic' = 'MIC', 
  'mic_tested' = 'Tested', 'comments' = 'X__1')

# FIX DRUG SYNONYMS ------------------------------------------------------------------
# ! was only run once - to update the respective files - should probably write a script for such 
# cases next time

# first problem: map the drugnames to drugname_typaslab
# View(tibble(mics$drugname, mics$drugname %in% drug_synonyms$drug_synonym))

# side note: all the information in 'drug_synonyms' is coming from 
# 'dbsetup/data/drug_synonyms_v2.RData' - but we don't know anymore how this was created
if (FALSE) {
  e <- new.env()
  dbsetup_drug_synonyms_file <- "~/PROJEKTE/dbsetup/data/drug_synonyms_v2.RData"
  load(dbsetup_drug_synonyms_file, envir = e)
  e$drug_synonyms_v2
  
  new_synonyms <- c('Bile acids' = 'BILE', 'Triton X-100' = 'TRITONX', 'CHIR- 090' = 'CHIR090', 
    'Epigallocatechin gallate (EGCG)' = 'EGCG', 'EPINEPHRINE' = 'EPINEPHRINE', 'Epinephrine' = 
      'EPINEPHRINE', 'NOREPINEPHRINE' = 'NOREPINEPHRINE', 'Norepinephrine' = 'NOREPINEPHRINE', 
    'EtOH' = 'ETHANOL', 'Nitrofurnatoin' = 'NITROFURANTOIN', 'Paraquat dichloride' = 'PARAQUAT', 
    'Phenazine methosulfate (PMS)' = 'PMS', 'CCCP (Carbonyl cyanide 3-chlorophenylhydrazone)' = 
      'CCCP', 'Triclosan/Irgasan' = 'TRICLOSAN', 'Irgasan' = 'TRICLOSAN', 'Hydrogen peroxide' = 
      'PEROXIDE')
  
  new_synonyms <- data.frame(drug_synonym = names(new_synonyms), drugname_typaslab = new_synonyms, 
    stringsAsFactors = FALSE)
  
  e$drug_synonyms_v2 <- rbind(e$drug_synonyms_v2, new_synonyms, stringsAsFactors = FALSE)
  # View(tibble(mics$drugname, mics$drugname %in% e$drug_synonyms_v2$drug_synonym))
  
  drug_synonyms <- e$drug_synonyms_v2
  
  # update tables and data:
  save(drug_synonyms, file = dbsetup_drug_synonyms_file)
  write_delim(drug_synonyms, delim = ";", path = drug_synonyms_file)
}

# CONTINUE ------------------------------------------------------------------------------------

# exclude cases that don't appear:
mics <- semi_join(mics, drug_synonyms, by = c('drugname_typaslab' = 'drug_synonym'))
mics$drugname_typaslab <- drug_synonyms$drugname_typaslab[match(x = mics$drugname_typaslab, 
  table = drug_synonyms$drug_synonym)]

mics$resistant <- ifelse(mics$mic == "resistant", TRUE, FALSE)
mics$mic_parsed <- parse_number(mics$mic, na = c("?", "-", "NA", "resistant"))
mics$mic_tested_parsed <- parse_number(mics$mic_tested, na = c("-", "NA", "yes"))
mics$mic_curated <- mics$mic_parsed
mics$mic_curated[!is.na(mics$mic_tested_parsed)] <- mics$mic_tested_parsed[!is.na(mics$mic_tested_parsed)]
mics <- select(mics, drugname_typaslab, mic, mic_parsed, mic_tested, mic_tested_parsed, 
  mic_curated, resistant, comments)
print(mics, n = 100)

corrections <- c("A22" = 2, "TRICLOSAN" = 0.8, "EPINEPHRINE" = NA, "NOREPINEPHRINE" = NA, 
  "NOVOBIOCIN" = 40, "CYCLOSERINED" = 10, "POLYMYXINB" = 0.5, "BACITRACIN" = NA, 
  "STREPTOMYCIN" = NA, "SDS" = 5, "BENZALKONIUM" = 20, "CCCP" = 50)

mics$mic_curated[match(names(corrections), mics$drugname_typaslab)] <- corrections
mics <- select(mics, drugname_typaslab, mic_curated, resistant, comments)

write_delim(mics, path = "./data/programmatic_output/MICs.csv", delim = ";")

