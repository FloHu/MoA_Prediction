## FROM COMMAND LINE -------------------------------
args <- commandArgs(trailingOnly = TRUE)
my_line <- as.numeric(args[1])

## PACKAGES -------------------------------
ipak <- function(pkg){
  suppressMessages(suppressWarnings({
    pkg = as.character(substitute(pkg))
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    invisible(sapply(pkg, require, character.only = TRUE))
  }))
}

ipak(MASS)
ipak(iterators)
ipak(foreach)
ipak(rpart)
ipak(xgboost)
ipak(glmnet)
ipak(tidyverse)
ipak(mlr)
ipak(parallelMap)

# parallel::detectCores()
parallelStartMulticore(cpus = 12)


## DATA AND PATHS -------------------------------
exportdir <- "/scratch/typas/cluster_run_results/mc_2019"
mc <- readRDS("mc.rds")


## FUNCTION DEFINITIONS -------------------------------
make_my_task <- function(dfm, targetvar = "process_broad", blockvar = NULL) {
  # dfm = drug-feature matrix
  # this function just removes drugname_typaslab, conc; then defines 
  # process_broad as target column for a task
  # can provide a blocking variable 
  dfm <- as.data.frame(dfm)
  if (!is.null(blockvar)) blocks <- make_blocks(dfm, blockvar = blockvar)
  
  to_remove <- c("drugname_typaslab", "conc")
  df <- dfm[, !colnames(dfm) %in% to_remove]
  
  le_task <- makeClassifTask(data = df, target = targetvar, 
    blocking = if (!is.null(blockvar)) blocks else NULL)
  
  # annoying thing about mlr task is that there is no option to save metadata 
  # about the task 
  # so we generate our own custom slot
  le_task$data_complete <- as_tibble(dfm)
  
  return(le_task)
}

make_blocks <- function(dfr, blockvar) {
  # takes a data frame where observations belonging together are next to each 
  # other in the column blockvar, then returns a factor containing the blocks 
  # (for mlr)
  rl <- rle(dfr[[blockvar]])
  blocks <- factor(rep(seq_along(rl$lengths), rl$lengths))
  return(blocks)
}

make_cvinst_blocked_stratified <- 
  function(task_data_all_cols, mlr_task, folds, strat_var) {
    # takes a task object from mlr and produces a CV instance with folds folds 
    # from it that is (i) stratified according to strat_var and (ii) blocked 
    # (blocking is already specified in the task)
    # task_data_all_cols: should contain also drugname_typaslab and conc to check if 
    # blocking worked correctly
    # mlr currently doesn't support this (specifying stratify = TRUE in 
    # makeResampleDesc() will throw an error in subsequent makeResampleInstance() 
    # functions if combined with blocked tasks)
    
    # ! intended for use with multiclass models, only ensures that training/test 
    # set are not empty
    stopifnot(is(mlr_task, "Task"))
    if (is.null(mlr_task$blocking)) {
      stop("Blocking has not been set for passed mlr task")
    }
    
    # usually we would set stratify = TRUE here
    rdesc <- mlr::makeResampleDesc("CV", iters = folds)
    data <- mlr::getTaskData(mlr_task)
    counter <- 0
    
    # make sure there is at least one observation in the test set of each class 
    # and in each train-test set split - this is our "stratification"
    try_again <- TRUE
    while (try_again) {
      # generate resampling instance
      rin <- mlr::makeResampleInstance(rdesc, mlr_task)
      # and check stratification
      for (split in seq_len(folds)) {
        tab <- table(data[[strat_var]][rin$test.inds[[split]]])
        # prevent empty classes in the test set
        if (any(tab == 0)) break
        # but also prevent empty classes in the training set
        tab <- table(data[[strat_var]][rin$train.inds[[split]]])
        if (any(tab == 0)) break
        try_again <- FALSE
      }
      counter <- counter + 1
    }
    
    # make sure that blocking worked correctly: drugs should not overlap 
    # between training and test set
    train_test_intersects <- 
      purrr::map2(rin$train.inds, rin$test.inds, function(train, test) {
        intersect(task_data_all_cols$drugname_typaslab[train], 
          task_data_all_cols$drugname_typaslab[test])
      })
    
    if (!all(sapply(train_test_intersects, length) == 0)) {
      stop("Something went wrong - drugs were found both in the train and 
        in the test set")
    }
    
    message("Generated blocked stratified resampling instance after ", counter, 
      " attempts")
    return(rin)
    }

fit_model_container_row <- function(row) {
  # takes a row from matrix container and runs repeated nested CV 
  # row = a row from matrix_container
  stopifnot({
    is(row, "data.frame")
    nrow(row) == 1
    !any(map_lgl(row, is.null))
    # name checking routine
  })
  
  row_params <- map(row[1, ], 1)
  n_rep <- length(row_params$resamp_instance)
  run_result <- list()
  
  # for each cv instance
  # for (cvrep in seq_len(n_rep)) {
  foreach(cvrep = icount(n_rep)) %do% {
    cat("This is nested CV, repetition", cvrep, "\n")
    
    cvinst <- row_params$resamp_instance[[cvrep]]
    train_test_pairs <- transpose(list(cvinst$train.inds, cvinst$test.inds))
    
    run_outer_fold <- list()
    # for each pair:
    # for (n in seq_len(length(train_test_pairs))) {
    foreach(n = icount(length(train_test_pairs))) %do% {
      cat("This is split ", n, " of CV repetition ", cvrep, "\n")
      
      fold_name <- paste0("outer_fold_", n)
      pair <- train_test_pairs[[n]]
      # make sure the pair is really a pair and that no observations 
      # occur in both training and test set
      stopifnot({
        length(pair) == 2
        length(pair[[1]]) > length(pair[[2]])
        n_distinct(flatten_dbl(pair)) == length(flatten_dbl(pair))
      })
      
      training_data <- 
        row_params$drug_feature_matrices[pair[[1]], ] %>%
        arrange(drugname_typaslab) # for correct blocking
      
      task_train_outer <- make_my_task(dfm = training_data, 
        blockvar = "drugname_typaslab", targetvar = "process_broad")
      
      test_data <- row_params$drug_feature_matrices[pair[[2]], ]
      
      task_test_outer <- make_my_task(dfm = test_data, 
        targetvar = "process_broad", blockvar = NULL)
      
      lrn_train_outer <- makeLearner(cl = row_params$fitted_model, predict.type = "prob")
      # check if this changes depending on algorithm
      # in principle the learner could also ship with the matrix_container ...
      lrn_train_outer_wrapped <- 
        makeFilterWrapper(learner = lrn_train_outer, fw.method = "variance", 
          fw.perc = 0.05)
      
      # tuning by cross-validation ("inner loop"), returns tuned learner, 
      # trained learner, and opt_path
      tune_res <- perform_tuning(.learner = lrn_train_outer_wrapped, 
        .data = training_data, .task = task_train_outer, .nfolds = 8, 
        .paramset = row_params$hyperparam_grid, .ctrl = makeTuneControlGrid())
      
      # predict test set
      prediction <- predict(tune_res$model, task = task_test_outer)
      stopifnot(all(prediction$data$truth == task_test_outer$data_complete$process_broad))
      prediction$data <- cbind(prediction$data, 
        task_test_outer$data_complete[, c("drugname_typaslab", "conc")])
      
      # save results
      run_outer_fold[[fold_name]] <- 
        tune_res[c("model", "opt_path", "opt_pars", "opt_pars_perf")]
      run_outer_fold[[fold_name]][["prediction"]] <- prediction
    }
    
    rep_name <- paste0("nested_cv_", cvrep)
    run_result[[rep_name]] <- run_outer_fold
  }
  
  return(run_result)
}

perform_tuning <- function(.learner, .data, .task, .nfolds, .paramset, .ctrl) {
  # input: learner with a dataset, number of folds for cv instance generation
  # then runs cross-validation-based parameter tuning and returns tuned learner 
  # plus the optimisation path
  # returns tuned learner + optimisation path
  
  .rin <- suppressMessages(
    make_cvinst_blocked_stratified(task_data_all_cols = .data, 
      mlr_task = .task, folds = .nfolds, strat_var = "process_broad")
  )
  
  # perform tuning
  tuning_result <- tuneParams(learner = .learner, task = .task, 
    resampling = .rin, measures = list(mmce, kappa), par.set = .paramset, 
    control = .ctrl)
  
  tuned_lrn <- setHyperPars(learner = .learner, par.vals = tuning_result$x)
  opt_path <- generateHyperParsEffectData(tuning_result)
  trained_model <- mlr::train(learner = tuned_lrn, task = .task)
  
  return(list(
    tuned_lrn = tuned_lrn, 
    opt_path = opt_path, 
    opt_pars = tuning_result$x, 
    opt_pars_perf = tuning_result$y, 
    model = trained_model))
}


## GET GOING -------------------------------

begin <- Sys.time()
result_name <- 
  purrr::reduce(mc[my_line, c("fitted_model", "drug_dosages", "chemical_feats")], 
    paste, sep = "_")

result <- fit_model_container_row(mc[my_line, ])

saveRDS(object = result, file = file.path(exportdir, paste0(result_name, ".rds")))

parallelStop()
end <- Sys.time()
cat(format(begin), "; ", format(end), "; ", format(end - begin), "; ", 
  result_name, "; ", my_line, "\n", file = "timings.txt", append = TRUE)






