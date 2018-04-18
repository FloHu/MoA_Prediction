get_measure_from_thresholds <- function(thresholds, perf_func, data, targetclass, prefix = "prob.") {
   # takes a vector of thresholds, a function for a performance metric, a prediction data frame
   # and a target class to get a vector of the chosen performance metric for targetclass vs. all
   # the other classes
   # thresholds: a vector of probability thresholds
   # perf_func: a function to calculate a performance measure from a confusion matrix
   # data: a prediction data frame
   # targetclass: the class for which the performance should be calculated
   # prefix: the prefix for targetclass used in data to apply the thresholds to
   if ( (!is.numeric(thresholds)) | any(c(range(thresholds) < 0, range(thresholds) > 1)) ) {
      stop("Provided vector thresholds is not numeric or thresholds don't lie between 0 and 1.")
   }
   if ( ! is.function(perf_func)) {
      stop("perf_func is not a function.")
   }

   probcol <- paste0(prefix, targetclass)
   if (! probcol %in% colnames(data)) {
      stop("Probabilities for the target class do not exist in the provided data.")
   }

   res <- sapply(thresholds, function(x) {
      tmp_dfr <- pred_df_to_onevsall_cm(pred_df = data, targetclass = targetclass,
                                        probability_col = probcol, thresh = x)
      perf_func(tmp_dfr)
   })
   return(res)
}
