knitr::opts_chunk$set(cache = FALSE, collapse = TRUE)
Sys.setlocale("LC_ALL", "en_IE.UTF-8")
source("./R/ipak.R")

ipak(magrittr)
ipak(tidyverse)
ipak(mlr)
ipak(parallelMap)
ipak(parallel)
ipak(plotmo)
ipak(reshape2)
ipak(gplots)
ipak(gridExtra)
ipak(plotly)
ipak(ComplexHeatmap)
ipak(circlize)
ipak(googledrive)
ipak(viridis)
ipak(mmpf)
ipak(ggpubr)
ipak(ggrepel)

ipak(extrafont)
loadfonts()

# custom functions
walk(list.files("./R", pattern = "*.R", full.names = T), source)
mics <- read_delim("./data/programmatic_output/MICs.csv", delim = ";")
mode_of_action <- read_csv("./data/programmatic_output/drug_moa_userfriendly.csv")

# default colours for mode of action, lookup tables for renaming
moa_cols <- c(cell_wall = "#1b9e77", dna = "#d95f02",
  membrane_stress = "#7570b3", protein_synthesis = "#e7298a", 
  unknown = "#636363")

moa_cols2 <- c(prob.cell_wall = "#1b9e77", prob.dna = "#d95f02", 
  prob.membrane_stress = "#7570b3", prob.protein_synthesis = "#e7298a")

moas <- names(moa_cols)
main_moas <- c("cell_wall", "dna", "membrane_stress", "protein_synthesis")

moa_repl_v0 <- c(cell_wall = "Cell Wall", membrane_stress = "Membrane Stress",
  protein_synthesis = "Protein Synthesis", dna = "DNA")

moa_repl <- c(cell_wall = "Cell Wall", membrane_stress = "Membrane\nStress",
  protein_synthesis = "Protein\nSynthesis", dna = "DNA")

moa_repl2 <- c(prob.cell_wall = "Cell Wall", prob.dna = "DNA", 
  prob.membrane_stress = "Membrane Stress", 
  prob.protein_synthesis = "Protein Synthesis")

moa_repl3 <- c(prob.cell_wall = "Cell Wall", prob.dna = "DNA", 
  prob.membrane_stress = "Membrane\nStress", 
  prob.protein_synthesis = "Protein\nSynthesis")

classifier_repl <- c(classif.randomForest = "Random Forests", 
  classif.glmnet = "Lasso", classif.xgboost = "Boosted Trees")

classifier_repl2 <- c(classif.randomForest = "Random\nForests", 
  classif.glmnet = "Lasso", classif.xgboost = "Boosted\nTrees")


global_fpt_order <- 
  c("MUTS", "STFR", "DGKA", "CYDB", "ISPA", "FADD", "FIMF", "DADX", "RBFA", 
    "NAGA", "PAL_and_4_more", "DDLB", "SLT", "MRCB", "YCFM", "YADB", "HUPA", 
    "HUPB", "RYFC", "ASMA", "LPCA_and_3_more", "RECG", "RECA", "LOLB", "RECC")

# plotting themes
theme_set(theme_bw())

comparison_theme <- theme_bw() + 
  theme(
  line = element_line(size = 0.1), # in mm
  text = element_text(size = 5) # in pts, acc. to doc - but this doesn't seem to be the case? 
)

theme(text = element_text(size = 12))

poster_theme <- theme(
  line = element_line(size = 0.5), 
  axis.title = element_text(size = 16, family = "Arial"),
  axis.text = element_text(size = 14, family = "Arial"), 
  strip.text = element_text(size = 14, family = "Arial"), 
  plot.title = element_text(size = 20, family = "Arial"), 
  plot.subtitle = element_text(size = 16, family = "Arial"), 
  legend.text = element_text(size = 14, family = "Arial"), 
  legend.title = element_text(size = 16, family = "Arial"), 
  panel.grid = element_blank()
)

paper_theme <- theme(
  line = element_line(size = 0.5), 
  axis.title = element_text(size = 8, family = "Arial"),
  axis.text = element_text(size = 7, family = "Arial"), 
  strip.text = element_text(size = 7, family = "Arial"), 
  plot.title = element_text(size = 8, family = "Arial"), 
  plot.subtitle = element_text(size = 7, family = "Arial"), 
  legend.text = element_text(size = 7, family = "Arial"), 
  legend.title = element_text(size = 8, family = "Arial"), 
  panel.grid = element_blank()
)
