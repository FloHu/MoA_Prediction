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
  
  # Approach: "difference to next highest probability": Class positive is the 
  # positive class and has been assigned a probability by the model. One can 
  # look at the highest probability of the other classes and take the highest 
  # one. The higher the distance to that probability the more confident we are 
  # that the current condition does indeed belong to the positive class. This 
  # distance can also be negative. 
  preproc <- group_by(preproc, drugname_typaslab, conc) %>% 
    mutate(next_prob = max(prob.med[prob_for != positive])) %>% 
    filter(prob_for == positive) %>% 
    ungroup()
  
  preproc$prob_delta <- preproc$prob.med - preproc$next_prob
  
  return(preproc)
}

get_responses_mcl <- function(dfr, positive, thresh) {
  # takes output of prep_roc
  dfr$positive <- positive
  dfr$thresh <- thresh
  negative <- paste0("not_", positive)
  dfr$response <- ifelse(dfr$prob_delta >= thresh, positive, negative)
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
  thresholds = seq(from = -1, to = 1, by = 0.01)) {
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
  
  # we can use geometry to calculate the AUC for bits of the ROC curve and then 
  # sum those areas up (similar to an integral)
  # programmatically, we do this by creating a pointlist, which contains the 
  # tpr and fpr of each point, and a pointlist_plusone, which contains the 
  # following points. Then we can map2() over these pairs of points to 
  # calculate the areas covered by each of points and at the end we sum 
  # everything up
  pointlist <- apply(threshs[, c("tpr", "fpr")], 1, as.list)
  pointlist_plusone <- dplyr::lead(pointlist)
  auc <- map2_dbl(pointlist, pointlist_plusone, function(x, y) {
    if (!is.list(y)) return(0)
    # not simplified mathematically for clarity
    area <- (y$fpr - x$fpr) * x$tpr + (y$fpr - x$fpr) * (y$tpr - x$tpr) / 2
  }) %>% 
    sum()
  
  p <- ggplot(threshs, aes(x = fpr, y = tpr, colour = positive)) + 
    geom_path(size = 0.75) + 
    geom_abline(slope = 1, size = 0.3) + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
    scale_colour_manual("Target process", values = moa_cols, 
      labels = moa_repl) + 
    paper_theme + 
    theme(panel.grid = element_blank(), legend.position = c(0.7, 0.2), 
      legend.background = element_blank()) + 
    labs(y = "FPR (1-specificity)", x = "TPR (recall)")
  
  p2 <- ggplot(threshs, aes(x = tpr, y = ppv, colour = positive)) + 
    geom_path(size = 0.75) + 
    # geom_abline(slope = 1, size = 0.3) + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
    scale_colour_manual("Target process", values = moa_cols, 
      labels = moa_repl) + 
    paper_theme + 
    theme(panel.grid = element_blank(), legend.position = c(0.7, 0.2), 
      legend.background = element_blank()) + 
    labs(x = "PPV (precision)", y = "TPR (recall)")
  
  print(p)
  invisible(list(auc_plot = p, prec_recall_plot = p2, auc = auc))
}
