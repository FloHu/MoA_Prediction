source("./R/ipak.R")
source("./setup.R")

ipak("glmnet")

# custom functions
walk(list.files("./R", pattern = "*.R", full.names = T), source)

load("./data/the_matrix_allDrugs.RData")
load("./data/the_matrix_newDrugs.RData")

load("./data/the_matrix_allDrugs_top10pct.RData.RData")
load("./data/the_matrix_newDrugs_top10pct.RData.RData")

load(file = "data/Rep_Nest_CV_instance_allDrugs.RData")
load(file = "data/Rep_Nest_CV_instance_newDrugs.RData")


lasso_hyp_param = makeParamSet(
    makeDiscreteParam("alpha", values = c(1)) #Lasso penality 
)
lasso_tuning = makeTuneControlGrid()



result_lasso_newDrugs = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs, model = "classif.glmnet", rep_instance = Rep_Nest_CV_instance_newDrugs,
                                                            run_hyp_param = lasso_hyp_param, run_tuning = lasso_tuning, tuning_measure = mmce)

result_lasso_allDrugs = rep_nested_CV_run_4models(data_matrix = the_matrix_allDrugs, model = "classif.glmnet", rep_instance = Rep_Nest_CV_instance_allDrugs,
                                                            run_hyp_param = lasso_hyp_param, run_tuning = lasso_tuning, tuning_measure = mmce)


lasso_hyp_param2 = makeParamSet(
    makeDiscreteParam("alpha", values = c(1)),  #Lasso penality 
    makeDiscreteParam("s", values = c(0.005,0.1,0.15,0.2,0.25))
)

result_lasso_newDrugs_tuneS = rep_nested_CV_run_4models(data_matrix = the_matrix_newDrugs, model = "classif.glmnet", rep_instance = Rep_Nest_CV_instance_newDrugs,
                                                  run_hyp_param = lasso_hyp_param2, run_tuning = lasso_tuning, tuning_measure = mmce)
