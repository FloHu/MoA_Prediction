# ==============================================================================
#
#                           CLUSTER RUN 
#							NEW DRUGS LABELLING          
#							XGBT - NEW HYPERPARAMETER GRID
#							TEST ADD LASSO AND SUBSAMPLIMG AND OTHER STUFF
#                           16/05/2018
#
# ==============================================================================

filepath_for_export = "/scratch/typas/cluster_run_results/"
filepath_for_import = "/home/dubois/MoA_Prediction/"



load(file = paste(filepath_for_import, "data/the_matrix_newDrugs.RData", sep=""))

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
    makeDiscreteParam("nrounds", values = c(200,500,1000)), # number of trees
    makeDiscreteParam("max_depth", values = c(1,3,5)), #depth 1 = stump, other depths enable interaction
    makeDiscreteParam("alpha", values = c(1)),#LASSO !
    makeDiscreteParam("eta", values = c(0.1,0.3)),
    makeDiscreteParam("subsample", values = c(0.5,0.8,1)), #individual subsample
    makeDiscreteParam("max_delta_step", values = c(0,5,10)) #Help estimate weight if class are umbalanced, 0 by default
 )
bt_tuning = makeTuneControlGrid()


result_XGBT_newDrugs_completeLasso = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
                                                            run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_newDrugs_completeLasso, file = paste0(filepath_for_export, "result_XGBT_newDrugs_completeLasso.RData"))


# ==============================================================================

#Hyperparameter tuning
bt_hyp_param = makeParamSet(
    makeDiscreteParam("nrounds", values = c(200,500,1000)), # number of trees
    makeDiscreteParam("max_depth", values = c(1,3,5)), #depth 1 = stump, other depths enable interaction
    makeDiscreteParam("alpha", values = c(0)),
    makeDiscreteParam("eta", values = c(0.1,0.3)),
    makeDiscreteParam("subsample", values = c(0.5,0.8,1)), #individual subsample
    makeDiscreteParam("colsample_bytree", values = c(0.25,0.5,0.75,1)),
    makeDiscreteParam("max_delta_step", values = c(0,5,10)) 
 )
bt_tuning = makeTuneControlGrid()


result_XGBT_newDrugs_completeFeatSample = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs, model = "classif.xgboost", rep_instance = Rep_Nest_CV_instance_newDrugs,
                                                            run_hyp_param = bt_hyp_param, run_tuning = bt_tuning, tuning_measure = mmce)
save(result_XGBT_newDrugs_completeFeatSample, file = paste0(filepath_for_export, "result_XGBT_newDrugs_completeFeatSample.RData"))


# ==============================================================================
