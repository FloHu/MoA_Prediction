fit_multiclass_rf <- function(data, 
                              targetvar, 
                              blockingvar = NULL, 
                              paramlist, 
                              cols_to_exclude = NULL, 
                              glue_back = TRUE) 
{
  # take a data data frame, blocked by blockingvar with optional 
  # cols_to_exclude (which are then glued back to the $pred$data slot if 
  # glue_back is TRUE)
  # to be added: a resampling argument
  data[[targetvar]] <- factor(data[[targetvar]])
  if (!is.null(blockingvar)) {
    my_rle <- rle(data[[blockingvar]])
    my_blocks <- factor(rep(seq_along(my_rle$lengths), my_rle$lengths))
  }
  
  if (!is.null(cols_to_exclude)) {
    taskdata <- dplyr::select(data, -one_of(cols_to_exclude))
  } else {
    taskdata <- data
  }
  
  task <- makeClassifTask(data = taskdata, target = targetvar, 
    blocking = if (!is.null(blockingvar)) my_blocks else NULL)
  lrnr <- makeLearner("classif.randomForest", predict.type = "prob", 
    par.vals = paramlist)
  rin <- make_rep_ncv(task_data_all_cols = data, mlr_task = task, reps = 10, 
    folds = 8, strat_var = targetvar)
  res <- resample(learner = lrnr, task = task, measures = mmce, resampling = 
      rin, models = TRUE, keep.pred = TRUE)
  
  if (glue_back) {
    res$pred$data <- dplyr::bind_cols(data[res$pred$data$id, cols_to_exclude], 
      res$pred$data)
  }
  
  return(res)
}


fit_multiclass_rf_loo <- function(data, 
  targetvar, 
  blockingvar = NULL, 
  paramlist, 
  cols_to_exclude = NULL, 
  glue_back = TRUE) 
{
  # take a data data frame, blocked by blockingvar with optional 
  # cols_to_exclude (which are then glued back to the $pred$data slot if 
  # glue_back is TRUE)
  # to be added: a resampling argument
  data[[targetvar]] <- factor(data[[targetvar]])
  if (!is.null(blockingvar)) {
    my_rle <- rle(data[[blockingvar]])
    my_blocks <- factor(rep(seq_along(my_rle$lengths), my_rle$lengths))
  }
  
  if (!is.null(cols_to_exclude)) {
    taskdata <- dplyr::select(data, -one_of(cols_to_exclude))
  } else {
    taskdata <- data
  }
  
  task <- makeClassifTask(data = taskdata, target = targetvar, 
    blocking = if (!is.null(blockingvar)) my_blocks else NULL)
  lrnr <- makeLearner("classif.randomForest", predict.type = "prob", 
    par.vals = paramlist)
  rdesc <- makeResampleDesc("LOO")
  
  res <- resample(learner = lrnr, task = task, measures = mmce, resampling = 
      rdesc, models = TRUE, keep.pred = TRUE)
  
  if (glue_back) {
    res$pred$data <- dplyr::bind_cols(data[res$pred$data$id, cols_to_exclude], 
      res$pred$data)
  }
  
  return(res)
}
