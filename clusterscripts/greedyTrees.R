# ==============================================================================
#
#                           CLUSTER RUN 
#							SUPER GREEDY TEST HYPERPARAMETER GRIDS           
#							COMPARE RF AND XGBT
#                            25/04/2018
#
# ==============================================================================

filepath_for_export = "/home/dubois/MoA_Prediction/run_results/"
filepath_for_import = "/home/dubois/MoA_Prediction/"



load(file = paste(filepath_for_import, "data/the_matrix_allDrugs_top10pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/Rep_Nest_CV_instance_allDrugs.RData", sep=""))


library(tidyverse)
library(mlr)
if(!require(rpart)){
    install.packages("rpart")
}
if(!require(randomForest)){
    install.packages("randomForest")
}
if(!require(xgboost)){
    install.packages("xgboost")
}
library(xgboost)
library(rpart)
library(randomForest)

walk(list.files(paste(filepath_for_import, "R/", sep = ""), pattern = "*.R", full.names = T), source)


library("parallelMap")
parallelStartMulticore(cpus = 24)



# ==============================================================================
#                                 RANDOM FOREST
# ==============================================================================

treespace <- c(100,200,300,500,1000)
rf_tuning = makeTuneControlGrid()

# ==============================================================================

data_set = the_matrix_allDrugs_top10pct
if("drugname_typaslab" %in% colnames(data_set)){
    data_set = select(data_set, -drugname_typaslab)
}

p = ncol(data_set)
mtry_space <- floor(c(p, p*(3/4), p/2, p/4, sqrt(p)))

rf_hyp_param <- makeParamSet(
    makeDiscreteParam("ntree", values = treespace), 
    makeDiscreteParam("mtry", values = mtry_space)
)

#result_greedy_RF_10pc_allDrugs_tuneMMCE = rep_nested_CV_run_4models(data_matrix = data_set, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning, tuning_measure = mmce)
#save(result_greedy_RF_10pc_allDrugs_tuneMMCE, file = paste0(filepath_for_export, "result_greedy_RF_10pc_allDrugs_tuneMMCE.RData"))



# ==============================================================================
#                                  BOOSTING TREE
# ==============================================================================

#Hyperparameter tuning
bt_hyp_param = makeParamSet(
    makeDiscreteParam("nrounds", values = c(200,500,1000)), # number of trees
    makeDiscreteParam("max_depth", values = c(1,2,3,5,10)), #depth 1 = stump, 2 leaves for 4 classes, makes no sense 
    makeDiscreteParam("eta", values = c(0.01,0.1,0.5))#Shrinkage to prevent overfitting
)
bt_tuning = makeTuneControlGrid()


result_greedy_XGBT_10pc_allDrugs_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_allDrugs_top10pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_greedy_XGBT_10pc_allDrugs_tuneMMCE, file = paste0(filepath_for_export, "result_greedy_XGBT_10pc_allDrugs_tuneMMCE.RData"))