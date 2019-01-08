prep_roc <- function(resampled_multiclass, positive) {
  # takes a resampled instance for multi-class prediction + a class that 
  # has been defined as positive and returns the predictions accordingly in 
  # the long format
  # intended to be used for making a ROC-curve
  stopifnot(is(resampled_multiclass, "ResampleResult"))
  stopifnot(positive %in% resampled_multiclass$pred$data$truth)
  stopifnot(length(positive) == 1)
  
  # a bit of renaming and data type changing
  preproc <- 
    resampled_multiclass %>%
    melt_pred_data(model_type = "multiclass") %>%
    mutate(truth = as.character(truth), 
      predicted_prob = str_replace(string = predicted_prob, 
        pattern = "prob\\.", replacement = "")) %>%
    dplyr::rename(prob_for = predicted_prob) %>%
    select(drugname_typaslab, conc, truth, prob_for, prob.med)
  
  # map all classes to either positive or not_positive
  negative <- paste0("not_", positive)
  preproc$truth <- ifelse(preproc$truth == positive, positive, negative)
  preproc[, c("truth", "prob_for")] <- 
    lapply(preproc[, c("truth", "prob_for")], function(col) {
      col <- ifelse(col == positive, positive, negative)
      return(col)
    })
  
  # keep only the probabilities for the positive class, the one for the 
  # negative is implicit
  preproc <- preproc[preproc$prob_for == positive, ]
  
  return(preproc)
}

get_responses_mcl <- function(dfr, positive, thresh) {
  # takes output of prep_roc
  dfr$positive <- positive
  dfr$thresh <- thresh
  negative <- paste0("not_", positive)
  
  dfr$response <- ifelse(dfr$prob.med >= thresh, positive, negative)
  dfr$tp <- (dfr$response == dfr$truth) & (dfr$truth == positive)
  dfr$fn <- (dfr$response != dfr$truth) & (dfr$truth == positive)
  dfr$fp <- (dfr$response != dfr$truth) & (dfr$truth != positive)
  dfr$tn <- (dfr$response == dfr$truth) & (dfr$truth != positive)
  stopifnot(sum(unlist(dfr[, c("tp", "fn", "fp", "tn")])) == nrow(dfr))
  return(dfr)
}

get_metrics_mcl <- function(dfr) {
  stopifnot(length(unique(dfr$thresh)) == 1)
  dfr <- group_by(dfr, thresh, positive) %>%
    summarise(tp = sum(tp), fn = sum(fn), tn = sum(tn), fp = sum(fp), 
      tpr = tp / (tp + fn), fpr = fp / (fp + tn), ppv = tp / (tp + fp))
  return(dfr)
}

get_thresh_vs_perf_mcl <- function(resampled_multiclass, positive, 
  thresholds = seq(from = 0, to = 1, by = 0.01)) {
  preproc <- prep_roc(resampled_multiclass, positive = positive)
  thresh_vs_perf <- map_dfr(thresholds, function(.thresh) {
    get_metrics_mcl(get_responses_mcl(preproc, positive, .thresh))
  })
  thresh_vs_perf <- arrange(thresh_vs_perf, desc(thresh))
  return(thresh_vs_perf)
}

plot_roc_mcl <- function(resampled_multiclass, positive) {
  get_thresh_vs_perf_mcl(resampled_multiclass, positive = positive) %>%
    ggplot(aes(x = fpr, y = tpr)) + 
    geom_path(aes(group = 1, colour = positive)) + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}