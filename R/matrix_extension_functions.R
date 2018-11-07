get_auc_info <- function(resultsobj, moas = c("cell_wall", "dna", "membrane_stress", 
  "protein_synthesis")) {
  # extract and merge all our prediction data from one of our results object
  # do this separately for each MoA
  # can be used as input to generateThreshVsPerfData()
  my_dfr <- expand.grid(cvrep = names(resultsobj), moa_modelled = moas, stringsAsFactors = FALSE)
  my_dfr <- as_tibble(my_dfr)
  
  my_dfr$predobj <-
    map2(my_dfr$cvrep, my_dfr$moa_modelled, function(.cvrep, .moa) {
      # we don't want a separate AUC for each outer fold but we need the structure of the
      # results object to concatenate the results of the other folds (this works correctly,
      # I checked)
      pred_obj_shell <- resultsobj[[.cvrep]][["Outer fold 1"]][[paste0("prediction_", .moa)]]
      pred_obj_shell$data <- map_dfr(resultsobj[[.cvrep]], function(.x) {
        return(.x[[paste0("prediction_", .moa)]]$data) # fetch the data slot for each fold of a specific .cvrep
      })
      return(pred_obj_shell)
    })
  
  my_dfr$auc <- map_dbl(my_dfr$predobj, ~ ROCR::performance(mlr::asROCRPrediction(.x), 
    measure = "auc")@y.values[[1]]) # accessory functions for ROCR?
  
  my_dfr$part_auc_01 <- map_dbl(my_dfr$predobj, ~ ROCR::performance(mlr::asROCRPrediction(.x), 
    measure = "auc", fpr.stop = 0.1)@y.values[[1]])
  
  # add 'ThreshVsPerfData' information as a separte column:
  my_dfr$thresh_vs_perf <- map(my_dfr$predobj, generateThreshVsPerfData, 
    measures = list(fpr, tpr, ppv))
  # add additional information:
  my_dfr$thresh_vs_perf <- imap(my_dfr$thresh_vs_perf, function(.data, .index) {
    .data$data$cvrep <- my_dfr$cvrep[.index]
    .data$data$moa_modelled <- my_dfr$moa_modelled[.index]
    return(.data)
  })
  
  return(my_dfr)
}


get_pred_data <- function(resultsobj, matrix_container_row, moas = c("cell_wall", "dna", 
  "membrane_stress", "protein_synthesis")) {
  # depends on functions prediction_merger and pred_extractor
  ## TO DO: proper error handling
  match.fun(prediction_merger)
  match.fun(pred_extractor)
  my_dfr <- map_dfr(moas, prediction_merger, resultsobj = resultsobj, extractorfunc = pred_extractor)
  my_dfr <- as_tibble(my_dfr)
  
  # add drug name and concentration information
  my_dfr$drugname_typaslab <- (matrix_container_row$drug_feature_matrices[[1]])$drugname_typaslab[my_dfr$id]
  my_dfr$conc <- (matrix_container_row$drug_feature_matrices[[1]])$conc[my_dfr$id]
  
  return(my_dfr)
}


add_info_to_pred_data <- function(pred_data, most_ia_mat) {
  # most_ia_mat = a drug_feature matrix containing the one concentration corresponding to most 
  # interactions
  most_ia_mat <- select(most_ia_mat, drugname_typaslab, conc, process_broad)
  
  joined <- select(pred_data, drugname_typaslab, conc) %>%
    left_join(most_ia_mat, by = c("drugname_typaslab" = "drugname_typaslab"), 
      suffix = c(".origin", ".mostias")) %>%
    mutate(conc.origin.is.mostias = (conc.origin == conc.mostias))
  
  pred_data$process_broad <- joined$process_broad
  pred_data$conc_mostias <- joined$conc.origin.is.mostias
  
  return(pred_data)
}


get_nsignif <- function(drugfeatmat, predDataTbl, the_matrix) {
  # function to add some more information to PredData: number of significant interactions, 
  # ranks of concentrations etc. 
  # drugfeatmat: from matrix_container_ext: the drug_feature matrix (used for modelling)
  # predDataTbl: from matrix_container_ext: the prediction table
  # the_matrix: from first notebook: the_matrix after the do(select_mutant(.)) step - contains 
  # sscores and qvalues for each drug-gene-concentration combination
  drugfeatmat_long <- 
    select(drugfeatmat, -process_broad) %>% 
    gather(3:ncol(.), key = "gene_synonym", value = "sscore")
  
  # by doing a join we can do a lookup which drug-concentration-gene combination is significant
  tab_nsignif <- left_join(drugfeatmat_long, the_matrix) %>%
    group_by(drugname_typaslab, conc) %>%
    summarise(n_signif = sum(significant, na.rm = TRUE))
  
  # and now join this back into the predDataTbl object - which is what we need for plotting
  # exclude drugs of category "other"
  tab_nsignif <- left_join(predDataTbl, tab_nsignif)
  tab_nsignif <- tab_nsignif[tab_nsignif$process_broad %in% 
      c("cell_wall", "dna", "membrane_stress", "protein_synthesis"), ]
  tab_nsignif$moa_modelled_is_truth <- tab_nsignif$moa_modelled == tab_nsignif$process_broad
  
  # in this case we're not interested in the probabilities for each repeat of the nested CV, only 
  # in the median
  tab_nsignif <- tab_nsignif %>%
    group_by(drugname_typaslab, conc, moa_modelled, moa_modelled_is_truth, n_signif, conc_mostias) %>%
    summarise(median_prob.moa = median(prob.moa))
  
  # get ranks of concentrations, of number of significant interactions, ...
  # but beware: conc with "highest number of significant interactions" was based on all features
  tab_nsignif <- 
    tab_nsignif %>%
    group_by(drugname_typaslab, moa_modelled) %>%
    mutate(conc_rank = rank(conc), 
      n_signif_rank = order(order(n_signif, conc)), # enables custom tie breaking as order(order(x)) == rank(x)
      n_conc = length(conc)) %>%
    ungroup(tmp)
  
  return(tab_nsignif)
}

# extractor functions get the predictions/ThreshVsPerfData from each test fold of each prediction
# object of each repat of the cross-validation
# the actual extraction is performed in prediction_merger()
# Measures as an argument with default values ? Thus same function for ROC or prec racal or mcc, etc VS threshold curves

# perf_extractor <- function(x, moa) {
# probably not needed any longer
#   # x = a prediction object
#   # moa = mode of action
#   perfs <- generateThreshVsPerfData(x[[paste0("prediction_", moa)]], measures = list(fpr, tpr, ppv))$data
#   return(perfs)
# }

# get_threshvsperf_data <- function(resultsobj, moas = c("cell_wall", "dna", "membrane_stress", 
#   "protein_synthesis")) {
#   # depends on prediction_merger and perf_extractor
#   ## TO DO: proper error handling
#   match.fun(prediction_merger)
#   match.fun(perf_extractor)
#   
#   map_dfr(moas, prediction_merger, resultsobj = resultsobj, extractorfunc = perf_extractor)
# }

pred_extractor <- function(x, moa) {
  # x = a prediction object
  # moa = mode of action
  x <- as.data.frame(x[[paste0("prediction_", moa)]])
  toreplace <- c(paste0("prob.", moa), paste0("prob.not_", moa))
  indstoreplace <- match(toreplace, names(x))
  names(x)[indstoreplace] <- c("prob.moa", "prob.not_moa")
  return(x)
}

prediction_merger <- function(resultsobj, moa, extractorfunc) {
  # resultsobj = run from a repeated nested CV
  # moa = one of our moas - check
  # extractorfunc = function (such as pred_extractor) that pulls out the stuff from the prediction object
  dfr_list <-
    map(resultsobj, function(.x) { # apply to each repeat of the nested CV
      map2(.x, names(.x), function(.x, .y) { # apply to each fold, record the name
        dfr <- extractorfunc(.x, moa)
        dfr$fold <- .y
        return(dfr)
      }) %>%
        bind_rows() # to rbind all the test folds
    })
  
  # collapse into a data frame
  dfr <- imap_dfr(dfr_list, function(.x, .y) {
    .x[["cvrep"]] <- .y
    return(.x)
  })
  
  # add moa information
  dfr[["moa_modelled"]] <- moa
  return(dfr)
}

make_filename <- function(matrix_container_row) {
  filename <- paste0(paste(unlist(matrix_container_row[, c("hyperparam_grid_name", "drug_dosages", 
    "feat_preselect", "chemical_feats")]), collapse = "_"), ".rds")
  return(filename)
}

extract_params_from_resultsobj <- 
  function(resobj, moa = c("dna", "cell_wall", "membrane_stress", "protein_synthesis"), param = c("ntree", "mtry")) {
    # takes a result object from matrix_container and extracts all the ntree/mtry values used 
    # across all folds and repeats
    moa <- match.arg(moa)
    param <- match.arg(param)
    
    get_param <- function(fold_contents, moa, param) {
      # each fold of each repeat contains model and prediction objects from which we 
      # can extract the corresponding hyperparameters
      getLearnerModel(fold_contents[[paste0("model_", moa)]])[["learner.model"]][[param]]
    }
    
    access_all_folds <- function(cvrepeat, FUN, ...) {
      # takes a cvrepeat and accesses all the folds using FUN
      # ... = additional arguments passed to FUN
      map(cvrepeat, FUN, ...)
    }
    
    # so what happens here?
    # the first call to map accesses all slots of the resobj, which are the nested CV repeats 
    # this accessing is done using access_all_folds + a function: access_all_folds will again 
    # call map so it will pass all the folds (hence the name) to FUN, which can then access 
    # the desired values
    params <- unlist(map(resobj, access_all_folds, FUN = get_param, param = param, moa = moa))
    return(params)
  }



