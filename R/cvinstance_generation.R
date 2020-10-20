get_outer_train_ids <- function(NCV_sampling, rep_nb){
  NCV_sampling[[rep_nb]]$outer$train.inds
}

get_inner_train_ids <- function(NCV_sampling, rep_nb, outer_nb){
  NCV_sampling[[rep_nb]]$inner[[outer_nb]]$train.inds
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
    rdesc <- mlr::makeResampleDesc("CV", iters = folds, blocking.cv = TRUE)
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


make_rep_ncv <- function(task_data_all_cols, 
                         mlr_task, 
                         reps, 
                         folds, 
                         strat_var, 
                         max_dev = 0.05) 
  {
  # make a resampling instance: not yet complete because it doesn't respect 
  # blocking and stratification: it's just to get an object of the correct 
  # structure
  # ! blocking needs to be specified in the task 
  rep_rin <- mlr::makeResampleInstance(mlr::makeResampleDesc(method = "RepCV", 
    reps = reps, folds = folds, blocking.cv = TRUE), mlr_task)
  
  # make the actual cv instances that are then used to overwrite the slots 
  # of rep_rin
  cv_insts <- replicate(
    reps, 
    make_cvinst_blocked_stratified(task_data_all_cols = task_data_all_cols, 
      mlr_task = mlr_task, folds = folds, strat_var = strat_var, 
      max_dev = max_dev), 
    simplify = FALSE
  )
  rep_rin$train.inds <- purrr::flatten(purrr::map(cv_insts, "train.inds"))
  rep_rin$test.inds <- purrr::flatten(purrr::map(cv_insts, "test.inds"))
  
  return(rep_rin)
  }

