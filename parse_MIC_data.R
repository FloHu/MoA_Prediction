library(readxl)
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

# MIC acriflavine? probably lower than 400, check papers, plus clear phenotype 
# already at 100
# CCCP: 25 after reviewing Ana's data
# mitomycin C: 1 after reviewing Ana's data
# ethidium bromide: changed to 150 based on doi: 10.1186/1754-1611-3-18
# nigericin: probably correct, see https://link.springer.com/chapter/10.1007/978-3-642-46051-7_45
# pyocyanin: changed to 20 based on Ana's data

corrections <- c("A22" = 2, "TRICLOSAN" = 0.8, "EPINEPHRINE" = NA, 
  "NOREPINEPHRINE" = NA, "NOVOBIOCIN" = 40, "CYCLOSERINED" = 10, 
  "POLYMYXINB" = 0.5, "BACITRACIN" = NA, "STREPTOMYCIN" = NA, "SDS" = 5, 
  "BENZALKONIUM" = 20, "CCCP" = 25, "FUSIDICACID" = 200, "STREPTOMYCIN" = 8, 
  "ACRIFLAVINE" = 100, "MITOMYCINC" = 1, "ETHIDIUMBROMIDE" = 150, 
  "PYOCYANIN" = 20)

mics$mic_curated[match(names(corrections), mics$drugname_typaslab)] <- corrections
mics <- select(mics, drugname_typaslab, mic_curated, resistant, comments)

write_delim(mics, path = "./data/programmatic_output/MICs.csv", delim = ";")

