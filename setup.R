knitr::opts_chunk$set(cache = F)
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

# custom functions
walk(list.files("./R", pattern = "*.R", full.names = T), source)

# default colours for mode of action
moa_cols <- c(cell_wall = "#1b9e77", dna = "#d95f02",
  membrane_stress = "#7570b3", protein_synthesis = "#e7298a")
moas <- names(moa_cols)
moa_repl <- c(cell_wall = "Cell Wall", membrane_stress = "Membrane Stress",
  protein_synthesis = "Protein Synthesis", dna = "DNA")

# plotting themes
theme_set(theme_bw())

comparison_theme <- theme(
  line = element_line(size = 0.1), # in mm
  text = element_text(size = 7),
  axis.text.x = element_text(angle = 45, hjust = 1)
)
