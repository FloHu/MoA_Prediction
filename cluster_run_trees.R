# ==============================================================================
#
#                   CLUSTER RUN FOR TREE BASED METHODS
#
# ==============================================================================

filepath_for_export = "/home/dubois/"
filepath_for_import = "/home/dubois/MoA_Prediction/"


load(file = paste(filepath_for_import, "the_matrix_hclust.RData", sep=""))
load(file = paste(filepath_for_import, "the_matrix_top10pct.RData", sep=""))
load(file = paste(filepath_for_import, "the_matrix_top5pct.RData", sep=""))

load(file = paste(filepath_for_import, "Rep_Nest_CV.RData", sep=""))


library(tidyverse)
library(mlr)
if(!require(adabag)){
    install.packages("adabag")
}
if(!require(rpart)){
    install.packages("rpart")
}
if(!require(randomForest)){
    install.packages("randomForest")
}
library(adabag)
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
    makeDiscreteParam("mfinal", values = c(500)), # or  seq(from = 100, to = 400, by = 20)
    makeDiscreteParam("maxdepth", values = c(3,4,5,6)), #depth 1 = stump, 2 leaves for 4 classes, makes no sense 
    makeDiscreteParam("coeflearn", values = "Zhu") #Best for multiclass problem
)
bt_tuning = makeTuneControlGrid()

# DONE !
#result_BT_10pc = rep_nested_CV_run(data_matrix = the_matrix_top10pct, model = "classif.boosting", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning)
#save(result_BT_10pc, file = paste(filepath_for_export, "result_BT_10pc.RData", sep = ""))
#result_BT_5pc = rep_nested_CV_run(data_matrix = the_matrix_top5pct, model = "classif.boosting", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning)
#save(result_BT_5pc, file = paste(filepath_for_export, "result_BT_5pc.RData", sep = ""))




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
result_RF_10pc = rep_nested_CV_run(data_matrix = the_matrix_top10pct, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning)
save(result_RF_10pc, file = paste(filepath_for_export, "result_RF_10pc.RData", sep = ""))

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
result_RF_5pc = rep_nested_CV_run(data_matrix = the_matrix_top5pct, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning)
save(result_RF_5pc, file = paste(filepath_for_export, "result_RF_5pc.RData", sep = ""))


# ============================================

data_set = the_matrix_hclust
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
result_RF_hclust = rep_nested_CV_run(data_matrix = the_matrix_hclust, model = "classif.randomForest", rep_instance = Rep_Nest_CV_instance, run_hyp_param = rf_hyp_param, run_tuning = rf_tuning)
save(result_RF_hclust, file = paste(filepath_for_export, "result_RF_hclust.RData", sep = ""))




#Run for BT with hclust, in the end because longer


result_BT_hclust = rep_nested_CV_run(data_matrix = the_matrix_hclust, model = "classif.boosting", rep_instance = Rep_Nest_CV_instance, run_hyp_param = bt_hyp_param, run_tuning = bt_tuning)
save(result_BT_hclust, file = paste(filepath_for_export, "result_BT_hclust.RData", sep = ""))
