fit_multiclass_rf <- function(learner, task, measures, rin, task_data_all_cols) {
  resampled_multiclass <- resample(learner = learner, task = task, 
    measures = measures, rin, models = TRUE, keep.pred = TRUE)
  
  resampled_multiclass$pred$data <- 
    dplyr::bind_cols(task_data_all_cols[resampled_multiclass$pred$data$id, c(1, 2)], 
      resampled_multiclass$pred$data) %>%
    dplyr::arrange(drugname_typaslab, conc, iter)
  
  return(resampled_multiclass)
}
