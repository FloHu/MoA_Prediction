# ==============================================================================
#
#                   CLUSTER RUN FOR XGBOOST AND NEURAL NETWORK
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

if(!require(xgboost)){
    install.packages("xgboost")
}
library(xgboost)

if(!require(nnet)){
    install.packages("nnet")
}
library(nnet)

walk(list.files(paste(filepath_for_import, "R/", sep = ""), pattern = "*.R", full.names = T), source)


library("parallelMap")
parallelStartMulticore(cpus = 24 )


# ==============================================================================
#                                  NOTE ON TUNING
# ==============================================================================
#
#	In this run, multiple tuning descision based on different measure minimization are tested
#	By default the hyperparameter set chosen is the one which minimizes mmce
#	will also test with mmcc multiclass.aunp and multiclass au1p (equivalent to AUC)
#	save every run as a seperate result so as to compare is the difference in final performances is significant



# ==============================================================================
#                                  XGBOOSTING
# ==============================================================================

xgb_hyp_param = makeParamSet(
    makeDiscreteParam("nrounds", values = c(500, 750, 1000)), # number of trees
    makeDiscreteParam("max_depth", values = c(4,5,6,10)), #depth 1 = stump, 2 leaves for 4 classes, makes no sense 
    makeDiscreteParam("eta", values = c(0.1,0.2,0.3)), #Shrinkage to prevent overfitting
    makeDiscreteParam("lambda", values = c(0.1,0.5,1)) # regularization to prevent overfitting
)
xgb_tuning = makeTuneControlGrid()

# with 10 % top variance features

result_xgb_10pc_mmce = rep_nested_CV_run(data_matrix = the_matrix_top10pct, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = mmce)
save(result_xgb_10pc_mmce, file = paste(filepath_for_export, "result_xgb_10pc_mmce.RData", sep = ""))

result_xgb_10pc_mmcc = rep_nested_CV_run(data_matrix = the_matrix_top10pct, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = multiclass_mcc)
save(result_xgb_10pc_mmcc, file = paste(filepath_for_export, "result_xgb_10pc_mmcc.RData", sep = ""))

result_xgb_10pc_au1p = rep_nested_CV_run(data_matrix = the_matrix_top10pct, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = multiclass.au1p)
save(result_xgb_10pc_au1p, file = paste(filepath_for_export, "result_xgb_10pc_au1p.RData", sep = ""))

result_xgb_10pc_aunp = rep_nested_CV_run(data_matrix = the_matrix_top10pct, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = multiclass.aunp)
save(result_xgb_10pc_aunp, file = paste(filepath_for_export, "result_xgb_10pc_aunp.RData", sep = ""))

# with 5 % top variance features

result_xgb_5pc_mmce = rep_nested_CV_run(data_matrix = the_matrix_top5pct, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = mmce)
save(result_xgb_5pc_mmce, file = paste(filepath_for_export, "result_xgb_5pc_mmce.RData", sep = ""))

result_xgb_5pc_mmcc = rep_nested_CV_run(data_matrix = the_matrix_top5pct, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = multiclass_mcc)
save(result_xgb_5pc_mmcc, file = paste(filepath_for_export, "result_xgb_5pc_mmcc.RData", sep = ""))

result_xgb_5pc_au1p = rep_nested_CV_run(data_matrix = the_matrix_top5pct, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = multiclass.au1p)
save(result_xgb_5pc_au1p, file = paste(filepath_for_export, "result_xgb_5pc_au1p.RData", sep = ""))

result_xgb_5pc_aunp = rep_nested_CV_run(data_matrix = the_matrix_top5pct, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = multiclass.aunp)
save(result_xgb_5pc_aunp, file = paste(filepath_for_export, "result_xgb_5pc_aunp.RData", sep = ""))


# with features selection based on hclust

result_xgb_hclust_mmce = rep_nested_CV_run(data_matrix = the_matrix_hclust, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = mmce)
save(result_xgb_hclust_mmce, file = paste(filepath_for_export, "result_xgb_hclust_mmce.RData", sep = ""))

result_xgb_hclust_mmcc = rep_nested_CV_run(data_matrix = the_matrix_hclust, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = multiclass_mcc)
save(result_xgb_hclust_mmcc, file = paste(filepath_for_export, "result_xgb_hclust_mmcc.RData", sep = ""))

result_xgb_hclust_au1p = rep_nested_CV_run(data_matrix = the_matrix_hclust, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = multiclass.au1p)
save(result_xgb_hclust_au1p, file = paste(filepath_for_export, "result_xgb_hclust_au1p.RData", sep = ""))

result_xgb_hclust_aunp = rep_nested_CV_run(data_matrix = the_matrix_hclust, model = "classif.xgboost", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = xgb_hyp_param, run_tuning = xgb_tuning, tuning_measure = multiclass.aunp)
save(result_xgb_hclust_aunp, file = paste(filepath_for_export, "result_xgb_hclust_aunp.RData", sep = ""))



# ==============================================================================
#                                  NEURAL NETWORK
# ==============================================================================



nn_hyp_param = makeParamSet(
    makeDiscreteParam("size", values = c(1, 3, 5, 10, 15, 20)), # number of hidden layers
    makeDiscreteParam("decay", values = c(0.001,0.005,0.01,0.05,0.1)),
    makeDiscreteParam("MaxNWts", values = 10**c(5,6)),
    makeDiscreteParam("maxit", values = c(500,750,1000)) #if converges stop before maxit
)
nn_tuning = makeTuneControlGrid()


# with 10 % top variance features

result_nn_10pc_mmce = rep_nested_CV_run(data_matrix = the_matrix_top10pct, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = mmce, saveFeatImportance = T)
save(result_nn_10pc_mmce, file = paste(filepath_for_export, "result_nn_10pc_mmce.RData", sep = ""))

result_nn_10pc_mmcc = rep_nested_CV_run(data_matrix = the_matrix_top10pct, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = multiclass_mcc, saveFeatImportance = T)
save(result_nn_10pc_mmcc, file = paste(filepath_for_export, "result_nn_10pc_mmcc.RData", sep = ""))

result_nn_10pc_au1p = rep_nested_CV_run(data_matrix = the_matrix_top10pct, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = multiclass.au1p, saveFeatImportance = T)
save(result_nn_10pc_au1p, file = paste(filepath_for_export, "result_nn_10pc_au1p.RData", sep = ""))

result_nn_10pc_aunp = rep_nested_CV_run(data_matrix = the_matrix_top10pct, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = multiclass.aunp, saveFeatImportance = T)
save(result_nn_10pc_aunp, file = paste(filepath_for_export, "result_nn_10pc_aunp.RData", sep = ""))

# with 5 % top variance features

result_nn_5pc_mmce = rep_nested_CV_run(data_matrix = the_matrix_top5pct, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = mmce, saveFeatImportance = T)
save(result_nn_5pc_mmce, file = paste(filepath_for_export, "result_nn_5pc_mmce.RData", sep = ""))

result_nn_5pc_mmcc = rep_nested_CV_run(data_matrix = the_matrix_top5pct, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = multiclass_mcc, saveFeatImportance = T)
save(result_nn_5pc_mmcc, file = paste(filepath_for_export, "result_nn_5pc_mmcc.RData", sep = ""))

result_nn_5pc_au1p = rep_nested_CV_run(data_matrix = the_matrix_top5pct, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = multiclass.au1p, saveFeatImportance = T)
save(result_nn_5pc_au1p, file = paste(filepath_for_export, "result_nn_5pc_au1p.RData", sep = ""))

result_nn_5pc_aunp = rep_nested_CV_run(data_matrix = the_matrix_top5pct, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = multiclass.aunp, saveFeatImportance = T)
save(result_nn_5pc_aunp, file = paste(filepath_for_export, "result_nn_5pc_aunp.RData", sep = ""))


# with features selection based on hclust

result_nn_hclust_mmce = rep_nested_CV_run(data_matrix = the_matrix_hclust, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = mmce, saveFeatImportance = T)
save(result_nn_hclust_mmce, file = paste(filepath_for_export, "result_nn_hclust_mmce.RData", sep = ""))

result_nn_hclust_mmcc = rep_nested_CV_run(data_matrix = the_matrix_hclust, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = multiclass_mcc, saveFeatImportance = T)
save(result_nn_hclust_mmcc, file = paste(filepath_for_export, "result_nn_hclust_mmcc.RData", sep = ""))

result_nn_hclust_au1p = rep_nested_CV_run(data_matrix = the_matrix_hclust, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = multiclass.au1p, saveFeatImportance = T)
save(result_nn_hclust_au1p, file = paste(filepath_for_export, "result_nn_hclust_au1p.RData", sep = ""))

result_nn_hclust_aunp = rep_nested_CV_run(data_matrix = the_matrix_hclust, model = "classif.nnet", 
					rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, tuning_measure = multiclass.aunp, saveFeatImportance = T)
save(result_nn_hclust_aunp, file = paste(filepath_for_export, "result_nn_hclust_aunp.RData", sep = ""))


