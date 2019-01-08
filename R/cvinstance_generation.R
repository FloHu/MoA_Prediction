get_outer_train_ids <- function(NCV_sampling, rep_nb){
  NCV_sampling[[rep_nb]]$outer$train.inds
}

get_inner_train_ids <- function(NCV_sampling, rep_nb, outer_nb){
  NCV_sampling[[rep_nb]]$inner[[outer_nb]]$train.inds
}

instance_creation <- function(dataset, printTest = F, nFoldsOuter = 8, nFoldsInner = 8){
  # takes a dataset = a drug feature matrix from matrix_container 
  # first, an outer CV instance is generated, which will later be used for performance evaluation 
  # second, for each training set of the outer CV instance another CV instance is generated, which 
  # will later be used for parameter tuning
  # returns a list with two components: (i) the outer CV instance, (ii) a list with one CV 
  # instance for each training set of the outer CV instance
  # stratification is defined in makeClassifTask() using "process_broad"
  dataset = select(dataset, -drugname_typaslab, -conc)
  dataset <- as.data.frame(dataset)
  sampling_method = makeResampleDesc(method = "CV", iters = nFoldsOuter, stratify = TRUE)
  
  # outer instance
  outer_task =  makeClassifTask(data = dataset, target = "process_broad")
  sample_instance_outer =  makeResampleInstance(sampling_method, outer_task)
  
  if (printTest) {
    # Testing stratification of outer instance
    for(i in 1:nFoldsOuter){
      cat("Testing stratification of outer instance.\nFold", i, "- Train and test sets\n" )
      ind = sample_instance_outer$train.inds[[i]]
      print(table(dataset[ind, "process_broad"]))
      ind = sample_instance_outer$test.inds[[i]]
      print(table(dataset[ind, "process_broad"]))
    }
    cat("\n\n\n#=================================================================\n\n\n")
  }
  
  # defining inner instances
  sampling_method = makeResampleDesc(method = "CV", iters = nFoldsInner, stratify = TRUE)
  sample_instance_inner = list()
  for(i in 1:nFoldsOuter){
    # New task using a subset of the whole dataset defined by outer training set ID
    inner_task = makeClassifTask(data = dataset[sample_instance_outer$train.inds[[i]], ], 
      target = "process_broad")
    sample_instance_inner[[i]] = makeResampleInstance(sampling_method, inner_task)
  }
  
  if(printTest){
    # Testing stratification of inner instances
    # NOTE: indices of individuals in the inner fold are the indices of the vector of indices of 
    # the corresponding outer fold 
    for (j in 1:nFoldsOuter ){
      outer_ind = sample_instance_outer$train.inds[[j]]
      cat(rep("#", 80), "\n", sep = "")
      cat("OUTER FOLD", j, ": training indices: \n")
      print(table(dataset[outer_ind, "process_broad"] ))
      cat("INNER FOLDS:\n")
      for (i in 1:nFoldsInner){
        train_inds <- sample_instance_inner[[j]]$train.inds[[i]]
        cat("Training set inner split", i, ":\n")
        print(table(dataset[outer_ind[train_inds], "process_broad"]))
        
        cat("Test set inner split", i, ":\n")
        test_inds <- sample_instance_inner[[j]]$test.inds[[i]]
        print(table(dataset[outer_ind[test_inds], "process_broad"]))
        cat(rep("-", 60), "\n", sep = "")
      }
    }
  }
  
  nested_CV_instance = list(outer = sample_instance_outer, inner = sample_instance_inner)
  return(nested_CV_instance)
}

RepNCV_instance_map_drugname <- function(instance_oneDrug, dataset_oneDosage, dataset_allDosage) {
  # This function takes a resampling instance made for one dosage (based on dataset_oneDosage) and 
  # from this generates a resampling instance for all dosages - based on dataset_allDosage - made 
  # by mapping the drugs with one dosage to the drugs with all dosages.
  # In this way, blocking is automatically respected.

  # the new resampling object will have the same structure as the original one
  Rep_Nest_CV_instance_allDosage = instance_oneDrug
  
  # for each nested cv instance repetition:
  for(i in 1:length(Rep_Nest_CV_instance_allDosage)) {
    # ====== OUTER CV MAPPING: ======
    # ======== Training sets ========
    # The indices of the training sets of the outer resampling are based on the drug_feature_matrix 
    # with one dosage = dataset_oneDosage
    # We want to map them to the indices of the dataset with all dosages.
    # To do so, we first retrieve the drug names that the indices of the training sets point to:
    outer_train_ids <- get_outer_train_ids(instance_oneDrug, rep_nb = i)
    drugNameSet_outer_train <- map(outer_train_ids, ~ dataset_oneDosage[.x, "drugname_typaslab"])
    
    # Then we take the drug_feature_matrix with all dosages (= dataset_allDosage) and use it to 
    # get the indices of all the drug names that the training set points to:
    allDosage_indices = map(drugNameSet_outer_train, function(x) {
      which(dataset_allDosage$drugname_typaslab %in% x$drugname_typaslab)
    })
    # and now replace the indices
    Rep_Nest_CV_instance_allDosage[[i]]$outer$train.inds = allDosage_indices 
    
    # ====== Test sets ======
    # Just use setdiff, what's not in the training set is in the test set
    allDosage_indices_test = map(allDosage_indices, function(x) {
      setdiff(seq_len(nrow(dataset_allDosage)), x)
    })
    # and again replace the indices
    Rep_Nest_CV_instance_allDosage[[i]]$outer$test.inds = allDosage_indices_test
    #Size mapping
    Rep_Nest_CV_instance_allDosage[[i]]$outer$size = nrow(dataset_allDosage)
    
    # INNER FOLDS RESAMPLING MAPPING
    for(outer in 1:length(outer_train_ids)) {
      # NOTE: NEVER MAP INNER IDS DIRECTLY TO THE DATASET, ALWAYS MAP THEM FIRST TO OUTER ID!
      # ====== Training sets ======
      # get all training indices of the inner CV for the corresponding CV repeat
      innerID = get_inner_train_ids(NCV_sampling = instance_oneDrug, rep_nb = i, outer_nb = outer)
      # but those innerID values are in fact indices of the following outer fold:
      outerID = get_outer_train_ids(NCV_sampling = instance_oneDrug, rep_nb = i)[[outer]]
      # so we need to map them to the outer IDs: the innerIDs subset the outerIDs and the resulting 
      # indices can then be actually used on the data
      innerID_to_rawData = map(.x = innerID, .f = function(x) outerID[x])
      
      # and now similar to above: first get the drug names that the training set refers to
      drugNameSet_inner_train = map(innerID_to_rawData, function(x) {
        dataset_oneDosage[x, "drugname_typaslab"]
      })
      
      # Then get the data set of all dosages that is still available to the inner CV by subsetting 
      # the table (dataset_allDosage) with the training indices of the correct fold:
      outer_IDs_alldosg <- get_outer_train_ids(NCV_sampling = Rep_Nest_CV_instance_allDosage, rep_nb = i)[[outer]]
      dataset_outer_train <- dataset_allDosage[outer_IDs_alldosg, "drugname_typaslab"]
      
      # same as above: by asking which drugs were actually used for training in each split of the 
      # inner CV we can map to all dosages:
      allDosage_indices = map(drugNameSet_inner_train, function(x) {
        which(dataset_outer_train$drugname_typaslab %in% x$drugname_typaslab)
      })
      
      # update indices
      Rep_Nest_CV_instance_allDosage[[i]]$inner[[outer]]$train.inds = allDosage_indices  
      
      # ====== Test sets ======
      allDosage_indices_test = map(allDosage_indices, function(x) {
        setdiff(seq_len(nrow(dataset_outer_train)), x)
      })
      Rep_Nest_CV_instance_allDosage[[i]]$inner[[outer]]$test.inds = allDosage_indices_test
      Rep_Nest_CV_instance_allDosage[[i]]$inner[[outer]]$size = nrow(dataset_outer_train)
    }
    
  }
  return(Rep_Nest_CV_instance_allDosage)
}

make_cvinst_blocked_stratified <- function(task_data_all_cols, 
                                           mlr_task, 
                                           folds, 
                                           strat_var, 
                                           max_dev = 0.05) {
  # takes a task object from mlr and produces a CV instance with folds folds 
  # from it that is (i) stratified according to strat_var and (ii) blocked 
  # according to block_var
  # task_data_all_cols: should contain also drugname_typaslab and conc to check if 
  # blocking worked correctly
  # mlr currently doesn't support this (specifying stratify = TRUE in 
  # makeResampleDesc() will throw an error in subsequent makeResampleInstance() 
  # functions if combined with blocked tasks)
  # Relative fractions in strat_var are allowed to differ by max_dev
  stopifnot(is(mlr_task, "Task"))
  if (is.null(mlr_task$blocking)) {
    stop("Blocking has not been set for passed mlr task")
  }
  
  rdesc <- mlr::makeResampleDesc("CV", iters = folds)
  data <- mlr::getTaskData(mlr_task)
  
  ref_fracs <- table(data[[strat_var]]) / length(data[[strat_var]])
  counter <- 0
  try_again <- TRUE
  
  while (try_again) {
    # generate resampling instance
    rin <- mlr::makeResampleInstance(rdesc, mlr_task)
    # and check stratification
    for (split in seq_len(folds)) {
      tab <- table(data[[strat_var]][rin$train.inds[[split]]])
      # prevent empty classes in the test set
      if (any(tab == 0)) break
      train_fracs <- tab / sum(tab)
      # make sure that relative fractions don't differ more than max_dev
      if (any(abs(ref_fracs - train_fracs) > max_dev)) break
      # if we get here it means we were successful
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
    reps = reps, folds = folds), mlr_task)
  
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
