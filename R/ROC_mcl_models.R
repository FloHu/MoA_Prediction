prep_roc <- function(pred_data, positive) {
  preproc <- 
    pred_data %>%
    melt_pred_data_mcl() %>%
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

get_thresh_vs_perf_mcl <- function(pred_data, positive, 
  thresholds = seq(from = 0, to = 1, by = 0.01)) {
  preproc <- prep_roc(pred_data, positive = positive)
  thresh_vs_perf <- map_dfr(thresholds, function(.thresh) {
    get_metrics_mcl(get_responses_mcl(preproc, positive, .thresh))
  })
  thresh_vs_perf <- arrange(thresh_vs_perf, desc(thresh))
  return(thresh_vs_perf)
}

plot_roc_mcl <- function(pred_data, positives) {
  # simulates a number of 1-vs-rest models for all the classes passed in 
  # positive by adding up all probabilities falling into the "not" category
  threshs <- map_dfr(positives, get_thresh_vs_perf_mcl, 
    pred_data = pred_data)
  ggplot(threshs, aes(x = fpr, y = tpr, colour = positive)) + 
    geom_path(size = 0.75) + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
    scale_colour_manual("Target process", values = moa_cols, 
      labels = names(moa_cols)) + 
    labs(x = "FPR (1-specificity)", y = "TPR (recall)")
}
