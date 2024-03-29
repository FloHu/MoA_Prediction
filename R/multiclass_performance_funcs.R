get_wide_confmat <- function(resample_result) {
  # input: a ResampleResult
  # output: a confusion matrix in long format for easy plotting
  if (!is(resample_result, "ResampleResult")) stop("`resample_result` is not a ResampleResult object")
  
  prd_obj <- resample_result$pred
  cm <- calculateConfusionMatrix(prd_obj)
  # we don't need the margins
  cm <- data.frame(cm$result[-nrow(cm$result), -ncol(cm$result)], 
    row.names = NULL, stringsAsFactors = FALSE)
  
  cm$true <- colnames(cm)
  cm <- gather(cm, -one_of("true"), key = "predicted", value = "n_obs")
  
  cm$predicted <- factor(cm$predicted)
  cm$true <- factor(cm$true, levels = rev(levels(cm$predicted)))
  
  cm$byclass_recall <- cm$n_obs / sapply(split(cm$n_obs, cm$true), sum)[cm$true]
  cm$byclass_recall <- cut(cm$byclass_recall, breaks = seq(from = 0, to = 1, by = 0.1))
  
  return(cm)
}

plot_wide_confmat <- function(
  wide_confmat, 
  toplot = c("recall", "precision"), 
  cols = RColorBrewer::brewer.pal(6, "Blues"), 
  save = FALSE, 
  file = NULL
  ) {
  toplot <- match.arg(toplot)
  
  p <- ggplot(wide_confmat, aes(x = predicted, y = true)) + 
    geom_tile(aes(fill = if (toplot == "recall") byclass_recall else byclass_prec)) + 
    geom_text(aes(label = n_obs), size = 2) + 
    labs(x = "Predicted class", y = "True class") + 
    scale_x_discrete(position = "top", labels = moa_repl) + 
    scale_y_discrete(labels = moa_repl) + 
    scale_fill_manual(toplot, values = cols, drop = FALSE) + 
    paper_theme + 
    theme(legend.key.size = unit(0.25, "cm"), panel.grid = element_blank())
  
  print(p)
  
  if (save) {
    if (is.null(file)) stop()
    ggsave(plot = p, file = file, width = 87, height = 60, units = "mm")
  }
  invisible(p)
}

plot_prob_calib <- function(pred_data, title, printplot = TRUE, save = FALSE, 
  file = NULL) {
  prob_calib <- melt_pred_data_mcl(pred_data)
  prob_calib$prob_bin <- cut(prob_calib$prob.med, breaks = seq(from = 0, to = 1, by = 0.1))
  prob_calib$prob_is_for <- str_extract(prob_calib$predicted_prob, 
    pattern = "cell_wall|dna|membrane_stress|protein_synthesis")
  prob_calib$truth <- as.character(prob_calib$truth)
  prob_calib <- 
    group_by(prob_calib, prob_is_for, prob_bin) %>%
    summarise(true_fraction = mean(truth == prob_is_for))
  prob_calib$prob_is_for <- recode(prob_calib$prob_is_for, !!!moa_repl)
  
  p <- ggplot(prob_calib, aes(x = as.numeric(prob_bin) - 0.05, y = true_fraction)) + 
    geom_point() + 
    geom_line(aes(group = prob_is_for)) + 
    geom_abline(slope = 1/10, linetype = "dotted") + 
    facet_wrap( ~ prob_is_for, ncol = 2) + 
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_x_continuous(breaks = seq(from = 0.5, to = 9.5, by = 1), 
      labels = as.character(10 * seq(from = 0, to = 9, by = 1))) + 
    labs(x = "Probability bin midpoint (10% bins)", y = "Observed event percentage", 
      title = "Probability calibration plot") + 
    comparison_theme
  
  if (printplot) suppressWarnings(print(p))
  
  if (save) {
    if (is.null(file)) stop("Must provide a filename")
    suppressWarnings(
      ggsave(filename = file, plot = p, width = 87, height = 90, units = "mm")
    )
  }
}
