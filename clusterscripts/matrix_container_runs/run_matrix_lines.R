# ==============================================================================
#
#                           CLUSTER RUN 
#							MATRIX CONTAINER       
#                           13/06/2018
#
# ==============================================================================

args = commandArgs(trailingOnly=TRUE)

startLine = args[1]
stopLine = args[2]

filepath_for_export = "/scratch/typas/cluster_run_results/matrix_container_result/"
filepath_for_import = "/home/dubois/MoA_Prediction/"



load(file = paste0(filepath_for_import, "data/matrix_container.RData"))

cat("******\n", startLine, "\n******\n")

library(tidyverse)
library(mlr)
if(!require(rpart)){
    install.packages("rpart")
}

if(!require(xgboost)){
    install.packages("xgboost")
}
if(!require(glmnet)){
    install.packages("glmnet")
}
library(xgboost)
library(rpart)


walk(list.files(paste(filepath_for_import, "R/", sep = ""), pattern = "*.R", full.names = T), source)


library("parallelMap")
parallelStartMulticore(cpus = 24)


for(line in startLine:stopLine){

	result = repeated_NCV_run_4models_container(data_container = matrix_container, line_number = line)

	result_name = paste(names(matrix_container[line,]$hyperparam_grid), 
						matrix_container[line,]$drug_dosages,
						matrix_container[line,]$feat_preselect,
						as.character(matrix_container[line,]$chemical_feats),
						sep = "_")
	saveRDS(object = result, file = paste0(filepath_for_export, result_name, ".rds"))	


}