# ==============================================================================
#
#                           CLUSTER RUN 
#           test Instance, dataset, feature slection and tuning effect on boosting tree
#
#                            24/04/2018
#
# ==============================================================================

filepath_for_export = "/home/dubois/MoA_Prediction/run_results/"
filepath_for_import = "/home/dubois/MoA_Prediction/"


load(file = paste0(filepath_for_import, "data/the_matrix.RData"))
load(file = paste(filepath_for_import, "data/the_matrix_top10pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_top5pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/Rep_Nest_CV_instance.RData", sep=""))

load(file = paste0(filepath_for_import, "data/the_matrix_allDrugs.RData"))
load(file = paste(filepath_for_import, "data/the_matrix_allDrugs_top10pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_allDrugs_top5pct.RData", sep=""))
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
#                                  BOOSTING TREE
# ==============================================================================

#Hyperparameter tuning
bt_hyp_param = makeParamSet(
    makeDiscreteParam("nrounds", values = c(200,500)), # number of trees
    makeDiscreteParam("max_depth", values = c(1,2,3)), #depth 1 = stump, 2 leaves for 4 classes, makes no sense 
    makeDiscreteParam("eta", values = c(0.01,0.1,0.5))#Shrinkage to prevent overfitting
)
bt_tuning = makeTuneControlGrid()





result_XGBT_5pc_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_top5pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning)
save(result_XGBT_5pc_tunePPV, file = paste0(filepath_for_export, "result_XGBT_5pc_tunePPV.RData"))

result_XGBT_5pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_top5pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_5pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_5pc_tuneMMCE.RData"))

result_XGBT_5pc_allDrugs_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_allDrugs_top5pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning)
save(result_XGBT_5pc_allDrugs_tunePPV, file = paste0(filepath_for_export, "result_XGBT_5pc_allDrugs_tunePPV.RData"))

result_XGBT_5pc_allDrugs_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_allDrugs_top5pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_5pc_allDrugs_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_5pc_allDrugs_tuneMMCE.RData"))



result_XGBT_10pc_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_top10pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning)
save(result_XGBT_10pc_tunePPV, file = paste0(filepath_for_export, "result_XGBT_10pc_tunePPV.RData"))

result_XGBT_10pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_top10pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_10pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_10pc_tuneMMCE.RData"))

result_XGBT_10pc_allDrugs_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_allDrugs_top10pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning)
save(result_XGBT_10pc_allDrugs_tunePPV, file = paste0(filepath_for_export, "result_XGBT_10pc_allDrugs_tunePPV.RData"))

result_XGBT_10pc_allDrugs_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_allDrugs_top10pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_10pc_allDrugs_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_10pc_allDrugs_tuneMMCE.RData"))




result_XGBT_wilcox_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, wilcoxSelection = T)
save(result_XGBT_wilcox_tunePPV, file = paste0(filepath_for_export, "result_XGBT_wilcox_tunePPV.RData"))

result_XGBT_wilcox_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce, wilcoxSelection = T)
save(result_XGBT_wilcox_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_wilcox_tuneMMCE.RData"))

result_XGBT_wilcox_allDrugs_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_allDrugs, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, wilcoxSelection = T)
save(result_XGBT_wilcox_allDrugs_tunePPV, file = paste0(filepath_for_export, "result_XGBT_wilcox_allDrugs_tunePPV.RData"))

result_XGBT_wilcox_allDrugs_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_allDrugs, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce, wilcoxSelection = T)
save(result_XGBT_wilcox_allDrugs_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_wilcox_allDrugs_tuneMMCE.RData"))



