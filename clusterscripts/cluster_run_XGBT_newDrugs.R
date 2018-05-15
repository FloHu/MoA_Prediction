# ==============================================================================
#
#                           CLUSTER RUN 
#							NEW DRUGS LABELLING          
#							XGBT - NEW HYPErPARAMETER GRID
#                           14/05/2018
#
# ==============================================================================

filepath_for_export = "/scratch/typas/cluster_run_results/"
filepath_for_import = "/home/dubois/MoA_Prediction/"



load(file = paste(filepath_for_import, "data/the_matrix_newDrugs.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_newDrugs_top5pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_newDrugs_top10pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_newDrugs_top15pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_newDrugs_top20pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_newDrugs_top25pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_newDrugs_top30pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_newDrugs_top40pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_newDrugs_top50pct.RData", sep=""))


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
    makeDiscreteParam("nrounds", values = c(200,500)), # number of trees
    makeDiscreteParam("max_depth", values = c(1,3,5)), #depth 1 = stump, other depths enable interaction
    makeDiscreteParam("lambda", values = c(0,0.01,0.1,1,10))#L2 Regularization (like a Ridge regression, 0 = no Ridge regression)
)
bt_tuning = makeTuneControlGrid()

# ==============================================================================


#
# 
# result_XGBT_newDrugs_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
# 								run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
# save(result_XGBT_newDrugs_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_newDrugs_tuneMMCE.RData"))
# 
# result_XGBT_newDrugs_tuneMCC = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
# 								run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mcc)
# save(result_XGBT_newDrugs_tuneMCC, file = paste0(filepath_for_export, "result_XGBT_newDrugs_tuneMCC.RData"))
# 
# 
# 
# result_XGBT_newDrugs_5pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top5pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
# 								run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
# save(result_XGBT_newDrugs_5pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_newDrugs_5pc_tuneMMCE.RData"))
# 
# result_XGBT_newDrugs_5pc_tuneMCC = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top5pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
# 								run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mcc)
# save(result_XGBT_newDrugs_5pc_tuneMCC, file = paste0(filepath_for_export, "result_XGBT_newDrugs_5pc_tuneMCC.RData"))
# 
# 
# 
# result_XGBT_newDrugs_10pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top10pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
# 								run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
# save(result_XGBT_newDrugs_10pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_newDrugs_10pc_tuneMMCE.RData"))
# 
# result_XGBT_newDrugs_10pc_tuneMCC = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top10pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
# 								run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mcc)
# save(result_XGBT_newDrugs_10pc_tuneMCC, file = paste0(filepath_for_export, "result_XGBT_newDrugs_10pc_tuneMCC.RData"))
# 
# 
# 
# result_XGBT_newDrugs_15pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top15pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
# 								run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
# save(result_XGBT_newDrugs_15pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_newDrugs_15pc_tuneMMCE.RData"))
# 
# result_XGBT_newDrugs_15pc_tuneMCC = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top15pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
# 								run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mcc)
# save(result_XGBT_newDrugs_15pc_tuneMCC, file = paste0(filepath_for_export, "result_XGBT_newDrugs_15pc_tuneMCC.RData"))
# 



# adding run on 15/05 only MMCE since the result are almost identical

# result_XGBT_newDrugs_20pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top20pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
#                                                               run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
# save(result_XGBT_newDrugs_20pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_newDrugs_20pc_tuneMMCE.RData"))


# result_XGBT_newDrugs_25pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top25pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
#                                                               run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
# save(result_XGBT_newDrugs_25pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_newDrugs_25pc_tuneMMCE.RData"))



result_XGBT_newDrugs_30pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top30pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
                                                               run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_newDrugs_30pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_newDrugs_30pc_tuneMMCE.RData"))


result_XGBT_newDrugs_40pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top40pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
                                                               run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_newDrugs_40pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_newDrugs_40pc_tuneMMCE.RData"))


result_XGBT_newDrugs_50pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs_top50pct, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
                                                               run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_newDrugs_50pc_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_newDrugs_50pc_tuneMMCE.RData"))

result_XGBT_newDrugs_all_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
                                                               run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_newDrugs_all_tuneMMCE, file = paste0(filepath_for_export, "result_XGBT_newDrugs_all_tuneMMCE.RData"))
