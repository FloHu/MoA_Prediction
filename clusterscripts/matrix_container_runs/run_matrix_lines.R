args = commandArgs(trailingOnly = TRUE)
startLine = as.numeric(args[1])
stopLine = as.numeric(args[2])

source("./ipak.R")
ipak(rpart)
ipak(xgboost)
ipak(glmnet)
ipak(tidyverse)
ipak(mlr)
ipak(foreach)
ipak(iterators)
ipak(parallelMap)

filepath_for_export <- "/scratch/typas/cluster_run_results/matrix_container_result"
matrix_container <- readRDS("matrix_container_2019.rds")

source("./repeated_NCV_run_4models_container.R")
source("./misc.R")
parallelStartMulticore(cpus = 24)

cat("******\nNow processing line ", startLine, "\n******\n")

for (line in startLine:stopLine) {
  result_name <- paste(names(matrix_container[line,]$hyperparam_grid), 
    matrix_container[line,]$drug_dosages,
    matrix_container[line,]$feat_preselect,
    as.character(matrix_container[line,]$chemical_feats),
    sep = "_")
  
  result <- repeated_NCV_run_4models_container(data_container = matrix_container, 
    line_number = line)

	saveRDS(object = result, file = file.path(filepath_for_export, 
	  paste0(result_name, ".rds")))
}

parallelStop()
