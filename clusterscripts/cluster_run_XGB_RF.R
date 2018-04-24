# ==============================================================================
#
#                   CLUSTER RUN FOR TREE BASED METHODS 
#
#                       run of 23/04/2018
#
# ==============================================================================

filepath_for_export = "/home/dubois/MoA_Prediction/run_results/"
filepath_for_import = "/home/dubois/MoA_Prediction/"


load(file = paste(filepath_for_import, "data/the_matrix_top10pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_top5pct.RData", sep=""))

load(file = paste(filepath_for_import, "data/Rep_Nest_CV.RData", sep=""))


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
parallelStartMulticore(cpus = 24 )



# ==============================================================================
#                                  BOOSTING TREE
# ==============================================================================

#Hyperparameter tuning
bt_hyp_param = makeParamSet(
    makeDiscreteParam("nrounds", values = c(600)), # number of trees
    makeDiscreteParam("max_depth", values = c(1,2,3,5,10)), #depth 1 = stump, 2 leaves for 4 classes, makes no sense 
    makeDiscreteParam("eta", values = c(0.01,0.1,0.3))#Shrinkage to prevent overfitting
)
bt_tuning = makeTuneControlGrid()



result_XGBT_10pc_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_top10pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning)
save(result_XGBT_10pc_tunePPV, file = paste0(filepath_for_export, "result_XGBT_10pc_tunePPV.RData"))
result_XGBT_5pc_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_top5pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning)
save(result_XGBT_5pc_tunePPV, file = paste0(filepath_for_export, "result_XGBT_5pc_tunePPV.RData"))


result_XGBT_10pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_top10pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_10pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_10pc_tuneMMCE.RData"))
result_XGBT_5pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_top5pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_5pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_5pc_tuneMMCE.RData"))



# ==============================================================================
#                              RANDOM FOREST
# ==============================================================================

treespace <- 200
rf_tuning = makeTuneControlGrid()

# ============================================

data_set = the_matrix_top10pct
if("drugname_typaslab" %in% colnames(data_set)){
    data_set = select(data_set, -drugname_typaslab)
}
if("process_broad" %in% colnames(data_set)){
    data_set = select(data_set, -process_broad)
}

p = ncol(data_set)
mtry_space <- floor(c(p, p*(3/4), p/2, p/4, sqrt(p)))

rf_hyp_param <- makeParamSet(
    makeDiscreteParam("ntree", values = treespace), 
    makeDiscreteParam("mtry", values = mtry_space)
)

result_RF_10pc_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_top10pct, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning)
save(result_RF_10pc_tunePPV, file = paste0(filepath_for_export, "result_RF_10pc_tunePPV.RData"))

result_RF_10pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_top10pct, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning, tuning_measure = mmce)
save(result_RF_10pc_tuneMMCE, file = paste0(filepath_for_export, "result_RF_10pc_tuneMMCE.RData"))

# ============================================

data_set = the_matrix_top5pct
if("drugname_typaslab" %in% colnames(data_set)){
    data_set = select(data_set, -drugname_typaslab)
}
if("process_broad" %in% colnames(data_set)){
    data_set = select(data_set, -process_broad)
}

p = ncol(data_set)
mtry_space <- floor(c(p, p*(3/4), p/2, p/4, sqrt(p)))

rf_hyp_param <- makeParamSet(
    makeDiscreteParam("ntree", values = treespace), 
    makeDiscreteParam("mtry", values = mtry_space)
)


result_RF_5pc_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_top5pct, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning)
save(result_RF_10pc_tunePPV, file = paste0(filepath_for_export, "result_RF_5pc_tunePPV.RData"))

result_RF_5pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_top5pct, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning, tuning_measure = mmce)
save(result_RF_10pc_tuneMMCE, file = paste0(filepath_for_export, "result_RF_5pc_tuneMMCE.RData"))
