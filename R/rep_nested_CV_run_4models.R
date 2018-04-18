#Script for repeated Nested Cross Validation run
#For a given dataset
#For a given model
#For a given HyperParamter grid
#For a given instance of resampling

rep_nested_CV_run_4models = function(data_matrix, model, rep_instance, run_hyp_param, 
                             run_tuning, saveFeatImportance = F, tuning_measure = ppv,
                             predict_type = "prob"){
    
    if(!require(foreach)){
        install.packages("foreach")
    }
    library(foreach)
    
    if("drugname_typaslab" %in% colnames(data_matrix)){
        data_matrix = select(data_matrix, -drugname_typaslab)
    }
    
    #Get number of repetition and CV folds
    n_rep = length(rep_instance)
    n_outer = length(rep_instance[[1]]$outer$train.inds)
    n_inner = length(rep_instance[[1]]$inner[[1]]$train.inds)
    

    #Final output object containing everything needed, models - prediction of all outer fold of all repetitions
    run_result = list()
    
    foreach(repetition_index = icount(n_rep)) %do% {
        #initialize outer folds output data for the repetition
        run_outer_fold = list()
        
        foreach(outerCV_ind = icount(n_outer)) %do% {
            
            outerCV_training_set = (get_outerCV(repetition_index, data = rep_instance))$train.ind[[outerCV_ind]]
            
            #define inner CV in relation to the outer fold used
            inner = get_innerCV(repetition_index, data = rep_instance)[[outerCV_ind]]
            
            #For or foreach ?
            #Seems better to put moa list inside the function, if arg would be less obvious when reading
            for (moa in c("dna", "cell_wall", "membrane_stress", "protein_synthesis")) {
                
                #create a nez matrix of all data with different process annotation
                custom_matrix = data_matrix %>% 
                                mutate(process_broad = replace(process_broad, process_broad != moa, paste0("not_", moa)))
                                
                #Thus the task should change because only a subset of the whole data should be used
                predictMoa = makeClassifTask(data = custom_matrix[ outerCV_training_set , ], target = "process_broad")
                
                #tuning hyperparameters based on the inner resampling
                resTuning = tuneParams(learner = model, task = predictMoa, resampling = inner, measures = tuning_measure,
                                       par.set = run_hyp_param, control = run_tuning)
                
                #use the best Hyperparams to create optimal learner
                run_learner = setHyperPars(makeLearner(cl = model, predict.type = predict_type), par.vals = resTuning$x)
            
                #We are now in the outer fold, all data should be used
                predictMoa = makeClassifTask(data = custom_matrix, target = "process_broad")
                
                #Model trained on a training subset
                model_outerCV = mlr::train(run_learner, predictMoa, subset = outerCV_training_set)
                
                outerCV_test_set = seq_along(1:nrow(custom_matrix))[-outerCV_training_set]
                pred_NCV = predict(model_outerCV, task = predictMoa, subset = outerCV_test_set)
                
                
                fold_name = paste0("Outer fold ", as.character(outerCV_ind))
                model_name = paste0("model_", moa)
                pred_name = paste0("prediction_", moa)
                run_outer_fold[[fold_name]][[model_name]] = model_outerCV 
                run_outer_fold[[fold_name]][[pred_name]] = pred_NCV
            }
            
            
            
            #if(saveFeatImportance){
            #    featimp = generateFeatureImportanceData(task = predictMoa, method = "permutation.importance",
            #                                            learner = run_learner)
            #    run_outer_fold[[fold_name]][["featuresImportance"]] = featimp 
            #}
            
        }
        
        rep_name = paste("Nested CV ", repetition_index, sep = "")
        run_result[[rep_name]] = run_outer_fold
    }
    

    return(run_result)
}
