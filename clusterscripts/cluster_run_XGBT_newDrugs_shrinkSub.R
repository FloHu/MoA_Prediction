# ==============================================================================
#
#                           CLUSTER RUN 
#							NEW DRUGS LABELLING          
#							XGBT - NEW HYPErPARAMETER GRID
#							TEST ADD SHRINKAGE AND SUBSAMPLIMG
#                           14/05/2018
#
# ==============================================================================

filepath_for_export = "/scratch/typas/cluster_run_results/"
filepath_for_import = "/home/dubois/MoA_Prediction/"



load(file = paste(filepath_for_import, "data/the_matrix_newDrugs.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_newDrugs_top10pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_newDrugs_top20pct.RData", sep=""))

load(file = paste(filepath_for_import, "data/Rep_Nest_CV_instance_newDrugs.RData", sep=""))


library(tidyverse)
library(mlr)
if(!require(rpart)){
    install.packages("rpart")
}

if(!require(xgboost)){
    install.packages("xgboost")
}
library(xgboost)
library(rpart)


walk(list.files(paste(filepath_for_import, "R/", sep = ""), pattern = "*.R", full.names = T), source)


library("parallelMap")
parallelStartMulticore(cpus = 24)

# ==============================================================================

#Hyperparameter tuning
bt_hyp_param = makeParamSet(
    makeDiscreteParam("nrounds", values = c(200,500,100)), # number of trees
    makeDiscreteParam("max_depth", values = c(1,3,5)), #depth 1 = stump, other depths enable interaction
    makeDiscreteParam("lambda", values = c(0,0.01,0.1,1,10)),#L2 Regularization (like a Ridge regression, 0 = no Ridge regression)
    makeDiscreteParam("eta", values = c(0.1,0.3)),
    makeDiscreteParam("subsample", values = c(0.5,0.8,1))
)
bt_tuning = makeTuneControlGrid()

# ==============================================================================



# result_XGBT_newDrugs_10pc_shrinkSub = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top10pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
#                                                              run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
# save(result_XGBT_newDrugs_10pc_shrinkSub, file = paste0(filepath_for_export, "result_XGBT_newDrugs_10pc_shrinkSub.RData"))



result_XGBT_newDrugs_20pc_shrinkSub = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top20pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
                                                            run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_newDrugs_20pc_shrinkSub, file = paste0(filepath_for_export, "result_XGBT_newDrugs_20pc_shrinkSub.RData"))


result_XGBT_newDrugs_all_shrinkSub = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
                                                            run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_newDrugs_all_shrinkSub, file = paste0(filepath_for_export, "result_XGBT_newDrugs_all_shrinkSub.RData"))

