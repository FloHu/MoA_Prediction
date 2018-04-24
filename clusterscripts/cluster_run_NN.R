# ==============================================================================
#
#                   CLUSTER RUN FOR NN METHODS 
#
#                       run of 20/04/2018
#
# ==============================================================================

filepath_for_export = "/home/dubois/MoA_Prediction/run_results/"
filepath_for_import = "/home/dubois/MoA_Prediction/"


load(file = paste(filepath_for_import, "data/the_matrix_top10pct.RData", sep=""))
load(file = paste(filepath_for_import, "data/the_matrix_top5pct.RData", sep=""))

load(file = paste(filepath_for_import, "data/Rep_Nest_CV.RData", sep=""))


library(tidyverse)
library(mlr)

if(!require(nnet)){
    install.packages("nnet")
}
library(nnet)

walk(list.files(paste(filepath_for_import, "R/", sep = ""), pattern = "*.R", full.names = T), source)


library("parallelMap")
parallelStartMulticore(cpus = 24 )



# ==============================================================================
#                                  NEURAL NETWORK
# ==============================================================================


nn_hyp_param = makeParamSet(
    makeDiscreteParam("size", values = c(1, 3, 5, 10)), # number of hidden layers
    makeDiscreteParam("decay", values = c(0.01,0.05,0.1)),
    makeDiscreteParam("MaxNWts", values = 10**c(7)),
    makeDiscreteParam("maxit", values = c(1000)) #if converges stop before maxit
)
nn_tuning = makeTuneControlGrid()


result_nn_10pc_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_top10pct, model = "classif.nnet", rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, saveFeatImportance = T)
save(result_nn_10pc_tunePPV, file = paste(filepath_for_export, "result_nn_10pc_tunePPV.RData", sep = ""))

result_nn_10pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_top10pct, model = "classif.nnet", rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, saveFeatImportance = T, tuning_measure = mmce)
save(result_nn_10pc_tuneMMCE, file = paste(filepath_for_export, "result_nn_10pc_tuneMMCE.RData", sep = ""))

result_nn_5pc_tunePPV = rep_nested_CV_run_4models(data_matrix = the_matrix_top5pct, model = "classif.nnet", rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, saveFeatImportance = T)
save(result_nn_5pc_tunePPV, file = paste(filepath_for_export, "result_nn_5pc_tunePPV.RData", sep = ""))

result_nn_5pc_tuneMMCE = rep_nested_CV_run_4models(data_matrix = the_matrix_top5pct, model = "classif.nnet", rep_instance = Rep_Nest_CV_instance, run_hyp_param = nn_hyp_param, run_tuning = nn_tuning, saveFeatImportance = T, tuning_measure = mmce)
save(result_nn_5pc_tuneMMCE, file = paste(filepath_for_export, "result_nn_5pc_tuneMMCE.RData", sep = ""))