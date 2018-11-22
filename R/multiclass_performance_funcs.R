get_wide_confmat <- function(resample_result) {
  # input: a ResampleResult
  # output: a confusion matrix in long format for easy plotting
  if (!is(resample_result, "ResampleResult")) stop("`resample_result` is not a ResampleResult object")
  supported_moas <- c("cell_wall", "dna", "membrane_stress", "protein_synthesis")

  prd_obj <- resample_result$pred
  cm <- calculateConfusionMatrix(prd_obj)
  cm <- data.frame(cm$result[-5, -5], row.names = NULL, stringsAsFactors = FALSE)
  
  if (!all(colnames(cm) %in% supported_moas)) stop("Only ", supported_moas, " are currently supported")
  
  cm$true <- colnames(cm)
  cm <- gather(cm, cell_wall:protein_synthesis, key = "predicted", value = "n_obs")
  
  cm$predicted <- factor(cm$predicted, levels = c("cell_wall", "dna", "membrane_stress", 
    "protein_synthesis"))
  cm$true <- factor(cm$true, levels = rev(levels(cm$predicted)))
  
  cm$byclass_recall <- cm$n_obs / sapply(split(cm$n_obs, cm$true), sum)[cm$true]
  cm$byclass_recall <- cut(cm$byclass_recall, breaks = seq(from = 0, to = 1, by = 0.1))
  
  return(cm)
}

plot_wide_confmat <- function(wide_confmat) {
  ggplot(wide_confmat, aes(x = predicted, y = true)) + 
    geom_tile(aes(fill = byclass_recall)) + 
    geom_text(aes(label = n_obs)) + 
    labs(x = "Predicted class, colours = recall/false-negative rate\n(normalisation by row)", 
      y = "True class") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0)) + 
    scale_x_discrete(position = "top") + 
    scale_fill_brewer()
}

plot_prob_calib <- function(resample_result) {
  # input: resample result
  if (!is(resample_result, "ResampleResult")) stop("`resample_result` is not a ResampleResult object")
  
  multicl_pred_melt <- melt_pred_data(resample_result, model_type = "multiclass")
  prob_calib <- multicl_pred_melt
  prob_calib$prob_bin <- cut(prob_calib$prob.med, breaks = seq(from = 0, to = 1, by = 0.1))
  levels(prob_calib$prob_bin)
  prob_calib$prob_is_for <- 
    str_extract(prob_calib$predicted_prob, pattern = "cell_wall|dna|membrane_stress|protein_synthesis")
  prob_calib$truth <- as.character(prob_calib$truth)
  prob_calib <- 
    group_by(prob_calib, prob_is_for, prob_bin) %>%
    summarise(true_fraction = mean(truth == prob_is_for))
  prob_calib
  
  p0 <- ggplot(prob_calib, aes(x = as.numeric(prob_bin) - 0.05, y = true_fraction)) + 
    geom_point() + 
    geom_line(aes(group = prob_is_for)) + 
    geom_abline(slope = 1/10, linetype = "dotted") + 
    facet_wrap( ~ prob_is_for, ncol = 2) + 
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_x_continuous(breaks = seq(from = 0.5, to = 9.5, by = 1), 
      labels = as.character(10 * seq(from = 0, to = 9, by = 1))) + 
    labs(x = "Probability bin midpoint (10% bins)", y = "Observed event percentage", 
      title = "Probability calibration plot for multi-class RF")
  p0
}


