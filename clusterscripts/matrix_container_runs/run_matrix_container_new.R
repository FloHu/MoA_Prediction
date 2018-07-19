# ==============================================================================
#
#                           CLUSTER RUN 
#							MATRIX CONTAINER       
#                           13/06/2018
#
# ==============================================================================


filepath_for_export = "/scratch/typas/cluster_run_results/matrix_container_result/"
filepath_for_import = "/home/dubois/MoA_Prediction/"



container = readRDS(file = paste0(filepath_for_import, "data/matrix_container_new.rds"))


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


for(line in 1:18){

	result = repeated_NCV_run_4models_container(data_container = container, line_number = line)

	result_name = paste(container[line, ]$drug_dosages,
						container[line, ]$fitted_model,
						container[line, ]$feat_preselect,
						as.character(container[line, ]$chemical_feats),
						sep = "_")
	saveRDS(object = result, file = paste0(filepath_for_export, result_name, ".rds"))	
}


