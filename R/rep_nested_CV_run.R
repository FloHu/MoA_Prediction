#Script for repeated Nested Cross Validation run
#For a given dataset
#For a given model
#For a given HyperParamter grid
#For a given instance of resampling

rep_nested_CV_run = function(data_matrix, model, rep_instance, run_hyp_param, run_tuning){
    
    if("drugname_typaslab" %in% colnames(data_matrix)){
        data_matrix = select(data_matrix, -drugname_typaslab)
    }
    
    #Get number of repetition and CV folds
    n_rep = length(rep_instance)
    n_outer = length(rep_instance[[1]]$outer$train.inds)
    n_inner = length(rep_instance[[1]]$inner[[1]]$train.inds)
    
    start_time = Sys.time()
    
    #Final output object containing everything needed, models - prediction of all outer fold of all repetitions
    run_result = list()

    for (repetition_index in 1:n_rep) {
        #initialize outer folds output data for the repetition
        run_outer_fold = list()
        
        for (outerCV_ind in 1:n_outer) {
            
            outerCV_training_set = (get_outerCV(repetition_index, data = rep_instance))$train.ind[[outerCV_ind]]
            
            #define inner CV in relation to the outer fold used
            inner = get_innerCV(repetition_index, data = rep_instance)[[outerCV_ind]]
            
            #Thus the task should change because only a subset of the whole data should be used
            predictMoa = makeClassifTask(data = data_matrix[ outerCV_training_set , ], target = "process_broad")
            
            #tuning hyperparameters based on the inner resampling
            resTuning = tuneParams(learner = model, task = predictMoa, resampling = inner,
                             par.set = run_hyp_param, control = run_tuning)
            
            #use the best Hyperparams to create optimal learner
            run_learner = setHyperPars(makeLearner(cl = model), par.vals = resTuning$x)
            
            #We are now in the outer fold, all data should be used
            predictMoa = makeClassifTask(data = data_matrix, target = "process_broad")
            
            #Model trained on a training subset
            model_outerCV = mlr::train(run_learner, predictMoa, subset = outerCV_training_set)
            
            outerCV_test_set = seq_along(1:nrow(data_matrix))[-outerCV_training_set]
            
            pred_NCV = predict(model_outerCV, task = predictMoa, subset = outerCV_test_set)
            
            #performance(pred_NCV, measures = list(mmce, kappa, multiclass_mcc))
            
            fold_name = paste("Outer fold ", as.character(outerCV_ind), sep = "")
            run_outer_fold[[fold_name]][["model"]] = model_outerCV 
            run_outer_fold[[fold_name]][["prediction"]] = pred_NCV
        }
        
        rep_name = paste("Nested CV ", repetition_index, sep = "")
        run_result[[rep_name]] = run_outer_fold
    }
    
    
    end_time = Sys.time()
    cat("Run finished - Time : ", end_time - start_time, " min", sep="")
    return(run_result)
}
