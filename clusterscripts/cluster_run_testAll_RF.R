# ==============================================================================
#
#                           CLUSTER RUN 
#           test Instance, dataset, feature slection and tuning effect on boosting tree
#
#                            25/04/2018
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

library(rpart)
library(randomForest)

walk(list.files(paste(filepath_for_import, "R/", sep = ""), pattern = "*.R", full.names = T), source)


library("parallelMap")
parallelStartMulticore(cpus = 24)



# ==============================================================================
#                                  BOOSTING TREE
# ==============================================================================

treespace <- c(200,500)
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

result_RF_10pc_allDrugs_tunePPV = rep_nested_CV_run_4models(data_matrix = data_set, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning)
save(result_RF_10pc_allDrugs_tunePPV, file = paste0(filepath_for_export, "result_RF_10pc_allDrugs_tunePPV.RData"))

result_RF_10pc_allDrugs_tuneMMCE = rep_nested_CV_run_4models(data_matrix = data_set, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning, tuning_measure = mmce)
save(result_RF_10pc_allDrugs_tuneMMCE, file = paste0(filepath_for_export, "result_RF_10pc_allDrugs_tuneMMCE.RData"))
 
# ==============================================================================

data_set = the_matrix_allDrugs_top5pct
if("drugname_typaslab" %in% colnames(data_set)){
    data_set = select(data_set, -drugname_typaslab)
}


p = ncol(data_set)
mtry_space <- floor(c(p, p*(3/4), p/2, p/4, sqrt(p)))

rf_hyp_param <- makeParamSet(
    makeDiscreteParam("ntree", values = treespace), 
    makeDiscreteParam("mtry", values = mtry_space)
)

result_RF_5pc_allDrugs_tunePPV = rep_nested_CV_run_4models(data_matrix = data_set, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning)
save(result_RF_5pc_allDrugs_tunePPV, file = paste0(filepath_for_export, "result_RF_5pc_allDrugs_tunePPV.RData"))

result_RF_5pc_allDrugs_tuneMMCE = rep_nested_CV_run_4models(data_matrix = data_set, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance_allDrugs, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning, tuning_measure = mmce)
save(result_RF_5pc_allDrugs_tuneMMCE, file = paste0(filepath_for_export, "result_RF_5pc_allDrugs_tuneMMCE.RData"))

# ==============================================================================

data_set = the_matrix_top5pct
if("drugname_typaslab" %in% colnames(data_set)){
    data_set = select(data_set, -drugname_typaslab)
}


p = ncol(data_set)
mtry_space <- floor(c(p, p*(3/4), p/2, p/4, sqrt(p)))

rf_hyp_param <- makeParamSet(
    makeDiscreteParam("ntree", values = treespace), 
    makeDiscreteParam("mtry", values = mtry_space)
)

result_RF_5pc_tunePPV = rep_nested_CV_run_4models(data_matrix = data_set, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning)
save(result_RF_5pc_tunePPV, file = paste0(filepath_for_export, "result_RF_5pc_tunePPV.RData"))

result_RF_5pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = data_set, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning, tuning_measure = mmce)
save(result_RF_5pc_tuneMMCE, file = paste0(filepath_for_export, "result_RF_5pc_tuneMMCE.RData"))


# ==============================================================================

data_set = the_matrix_top10pct
if("drugname_typaslab" %in% colnames(data_set)){
    data_set = select(data_set, -drugname_typaslab)
}


p = ncol(data_set)
mtry_space <- floor(c(p, p*(3/4), p/2, p/4, sqrt(p)))

rf_hyp_param <- makeParamSet(
    makeDiscreteParam("ntree", values = treespace), 
    makeDiscreteParam("mtry", values = mtry_space)
)



result_RF_10pc_tunePPV = rep_nested_CV_run_4models(data_matrix = data_set, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning)
save(result_RF_10pc_tunePPV, file = paste0(filepath_for_export, "result_RF_10pc_tunePPV.RData"))

result_RF_10pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = data_set, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning, tuning_measure = mmce)
save(result_RF_10pc_tuneMMCE, file = paste0(filepath_for_export, "result_RF_10pc_tuneMMCE.RData"))


