# ==============================================================================
#
#                           CLUSTER RUN 
#							MATRIX CONTAINER       
#                           BIAS INVESTIGATION
#
# ==============================================================================

filepath_for_export = "/scratch/typas/cluster_run_results/matrix_container_result/"
filepath_for_import = "/home/dubois/MoA_Prediction/"

load(file = paste0(filepath_for_import, "data/matrix_container.RData"))

library(tidyverse)
library(mlr)
if(!require(rpart)){
    install.packages("rpart")
}

if(!require(glmnet)){
    install.packages("glmnet")
}
library(rpart)


walk(list.files(paste(filepath_for_import, "R/", sep = ""), pattern = "*.R", full.names = T), source)


library("parallelMap")
parallelStartMulticore(cpus = 24)


# model Lasso_all(dosages)_all(features)_TRUE(with chemFeat)  
ref_model = matrix_container[146, ]
ref_model$tuning_measure[[1]] = ref_model$tuning_measure[[1]][c(2,3,1)]

ref_model_back = ref_model

# ========================= ELASTIC NET - VARIATION OF ALPHA =========================================



#ref_model$hyperparam_grid[[1]] = makeParamSet( makeDiscreteParam("alpha", values = c(0,0.2,0.5,0.7,1)), 
#                                           makeDiscreteParam("s", values = seq(from = 0.01, to = 10, length.out = 100)))

#res_elastic_net = repeated_NCV_run_4models_container(data_container = ref_model, line_number = 1)

#saveRDS(object = res_elastic_net, file = paste0(filepath_for_export, "result_elastic_net.rds"))	



# ========================= ONLY 4 MAIN MOA =========================================

# ref_model = matrix_container[164, ]


#ref_model$drug_feature_matrices[[1]] = ref_model$drug_feature_matrices[[1]] %>% filter(process_broad %in% c("dna", "cell_wall", "membrane_stress", "protein_synthesis"))

#RNCV_instance = list()
#dataset = ref_model$drug_feature_matrices[[1]]
#for(i in 1:10){
#    name = paste("NCV_", as.character(i), sep="")
#    RNCV_instance[[name]] = instance_creation(dataset = dataset, printTest = F)
#}

#ref_model$resamp_instance[[1]] = RNCV_instance


#res_main_moa = repeated_NCV_run_4models_container(data_container = ref_model, line_number = 1)

#saveRDS(object = res_main_moa, file = paste0(filepath_for_export, "res_main_moa.rds"))	


# ========================= EN WITH 1 DOSAGE TO COMPARE =========================================


ref_model = matrix_container[164, ]

ref_model$hyperparam_grid[[1]] = makeParamSet( makeDiscreteParam("alpha", values = c(0,0.2,0.5,0.7,1)), 
                                               makeDiscreteParam("s", values = seq(from = 0.01, to = 10, length.out = 100)))

res_EN_1dose = repeated_NCV_run_4models_container(data_container = ref_model, line_number = 1)

saveRDS(object = res_EN_1dose, file = paste0(filepath_for_export, "res_EN_1dose.rds"))	





# ========================= BOTH =========================================
# 
# ref_model = matrix_container[164, ]
# 
# ref_model$hyperparam_grid[[1]] = makeParamSet( makeDiscreteParam("alpha", values = c(0,0.2,0.5,0.7,1)), 
#                                          makeDiscreteParam("s", values = seq(from = 0.01, to = 10, length.out = 100)))
# 
# ref_model$drug_feature_matrices[[1]] = ref_model$drug_feature_matrices[[1]] %>% filter(process_broad %in% c("dna", "cell_wall", "membrane_stress", "protein_synthesis"))
# 
# 
# RNCV_instance = list()
# dataset = ref_model$drug_feature_matrices[[1]]
# for(i in 1:10){
#     name = paste("NCV_", as.character(i), sep="")
#     RNCV_instance[[name]] = instance_creation(dataset = dataset, printTest = F)
# }
# 
# ref_model$resamp_instance[[1]] = RNCV_instance
# res_main_moa_EN = repeated_NCV_run_4models_container(data_container = ref_model, line_number = 1)
# 
# saveRDS(object = res_main_moa_EN, file = paste0(filepath_for_export, "res_main_moa_EN.rds"))	
# 
# 
