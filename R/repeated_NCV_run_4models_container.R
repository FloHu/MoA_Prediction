repeated_NCV_run_4models_container = function(data_container, line_number, predict_type = "prob") {
  # Script for repeated Nested Cross Validation run
  # For a given dataset
  # For a given model
  # For a given HyperParameter grid
  # For a given instance of resampling
  # based on the respective fields of matrix_container
  if (line_number > nrow(data_container)) stop("Incorrect line number")

  data_matrix = data_container$drug_feature_matrices[[line_number]]
  model = data_container$fitted_model[[line_number]]
  rep_instance = data_container$resamp_instance[[line_number]]
  run_hyp_param = data_container$hyperparam_grid[[line_number]]
  tuning_measure = data_container$tuning_measure[[line_number]]
  
  data_matrix <- suppressWarnings(dplyr::select(data_matrix, -one_of(c("drugname_typaslab", "conc"))))

  # if (model == "classif.randomForest") {
  #     oobperf = makeMeasure("oobperf", minimize = TRUE, properties = c("classif", "classif.multi"), 
  #       fun = function(task, model, pred, feats, extra.args) {
  #         err = model$learner.model$err.rate
  #         err[nrow(err), "OOB"]
  #       })
  # }

  # Get number of repetition and CV folds
  n_rep = length(rep_instance)
  n_outer = length(rep_instance[[1]]$outer$train.inds)
  
  # Final output object containing everything needed
  run_result = list()
  
  # for each repetition of the nested CV:
  foreach(repetition_index = icount(n_rep)) %do% {
      # initialize outer folds output data for the repetition
      cat("REPETITION ", repetition_index, "\n")
      run_outer_fold = list()
      
      # for each train-test set split of each repetition:
      foreach(outerCV_ind = icount(n_outer)) %do% {
          cat("OUTER ", outerCV_ind, "\n")
          outerCV_training_set <- 
            get_outerCV(repeat_ind = repetition_index, data = rep_instance)$train.ind[[outerCV_ind]]
          
          # define inner CV in relation to the outer fold used
          inner = get_innerCV(repetition_index, data = rep_instance)[[outerCV_ind]]
          
          # fit a separate model for each MoA:
          for (moa in c("dna", "cell_wall", "membrane_stress", "protein_synthesis")) {
              cat("Mode of action:", moa, "\n")
              # create a new matrix of all data with different process annotation
              custom_matrix <- data_matrix
              custom_matrix$process_broad <- 
                ifelse(custom_matrix$process_broad %in% moa, custom_matrix$process_broad, paste0("not_", moa))
              custom_matrix$process_broad <- 
                factor(custom_matrix$process_broad, levels = c(moa, paste0("not_", moa)))
              custom_matrix = as.data.frame(custom_matrix)
              
              # Thus the task should change because only a subset of the whole data should be used
              predictMoa = makeClassifTask(data = custom_matrix[outerCV_training_set, ], target = "process_broad")
              
              # Use OOB error for tuning if model is random forest:
              if (FALSE) {#model == "classif.randomForest") {
                  # REMOVED because this doesn't respect blocking
                  # Tune wrapper is using an holdout resampling with a 0.99 split (1 observation as 
                  # test set, all others as training)
                  # It can be understood as "nearly no resampling"
                  # Measure evaluated is the OOB error extracted from the trained model object          
                  # tuned_lrn = makeTuneWrapper(
                  #   learner = makeLearner(cl = "classif.randomForest", predict.type = predict_type), 
                  #   par.set = run_hyp_param, control = makeTuneControlGrid(), measures = oobperf, 
                  #   resampling = makeResampleDesc("Holdout", split = 0.99)
                  # )
              } else {
                  # tuning hyperparameters based on the inner resampling
                  resTuning = tuneParams(learner = makeLearner(cl = model, predict.type = predict_type), 
                    task = predictMoa, resampling = inner, measures = tuning_measure, 
                    par.set = run_hyp_param, control = makeTuneControlGrid())
                  # use the best Hyperparams to create optimal learner.model
                  tuned_lrn = setHyperPars(makeLearner(cl = model, predict.type = predict_type), 
                    par.vals = resTuning$x)
              }
              
              # We are now in the outer fold, all data should be used
              predictMoa = makeClassifTask(data = custom_matrix, target = "process_broad")
              
              # Model trained on a training subset
              model_outerCV = mlr::train(tuned_lrn, predictMoa, subset = outerCV_training_set)
              
              outerCV_test_set = (1:nrow(custom_matrix))[-outerCV_training_set]
              pred_NCV = predict(model_outerCV, task = predictMoa, subset = outerCV_test_set)
              
              fold_name = paste0("Outer fold ", as.character(outerCV_ind))
              model_name = paste0("model_", moa)
              pred_name = paste0("prediction_", moa)
              run_outer_fold[[fold_name]][[model_name]] = model_outerCV 
              run_outer_fold[[fold_name]][[pred_name]] = pred_NCV
              # end of fitting a separate model for each MoA
          }
          # end of repeating the process for each training-test set split
      }
      rep_name = paste("Nested CV ", repetition_index, sep = "")
      run_result[[rep_name]] = run_outer_fold
      # all repetitions of the nested CV are done
  }
  return(run_result)
}
