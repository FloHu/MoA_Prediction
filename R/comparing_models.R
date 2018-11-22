melt_pred_data <- function(resample_result, model_type = c("onevsrest", "multiclass")) {
  # input: either a resample result or the pred_data data frame from matrix_container_ext
  model_type <- match.arg(model_type)
  
  if (model_type == "multiclass") {
    if (!is(resample_result, "ResampleResult")) stop("`resample_result` is not a ResampleResult object")
    pred_data <- resample_result$pred$data
  } else {
    pred_data <- resample_result
  }
  
  my_cols <- colnames(pred_data)
  
  # sanity check
  if (model_type == "onevsrest") {
    if(!(all(c("prob.moa", "prob.not_moa") %in% my_cols))) {
      stop("Passed data frame doesn't seem to fit the specified model type")
    }
  }

  # always needed
  if (!("conc" %in% my_cols)) {
    stop("Specified data frame doesn't contain a drug concentration")
  }
  
  # melt and average probabilities
  if (model_type == "multiclass") {
    tryCatch({
      melted <- tidyr::gather(pred_data, prob.cell_wall:prob.protein_synthesis, 
        key = "predicted_prob", value = "probability") %>%
        group_by(conc, predicted_prob, truth, drugname_typaslab) %>%
        summarise(prob.min = boxplot.stats(probability)$stats[2], 
          prob.max = boxplot.stats(probability)$stats[4], 
          prob.med = median(probability))
    }, 
      error = function(cnd) {
        stop("It looks like the data frame doesn't contain the four main modes of action.")
      })
  }
  
  if (model_type == "onevsrest") {
    tryCatch({
      pred_data$moa_modelled <- paste0("prob.", pred_data$moa_modelled)
      pred_data$truth <- pred_data$process_broad
      pred_data <- dplyr::rename(pred_data, probability = prob.moa, predicted_prob = moa_modelled)
      melted <- 
        pred_data %>%
        group_by(conc, predicted_prob, truth, drugname_typaslab) %>%
        summarise(prob.min = boxplot.stats(probability)$stats[2], 
          prob.max = boxplot.stats(probability)$stats[4], 
          prob.med = median(probability))
    }, 
      error = function(cond) {
        stop("There was a problem processing the one-vs-rest data frame")
      }
    )
  }
  
  return(melted)
}

compare_probabilities <- function(a, b, title = "") {
  # input: melted prediction data frames from melt_pred_data
  xlab <- deparse(substitute(a))
  ylab <- deparse(substitute(b))
  
  joined <- inner_join(a, b, by = c("conc", "predicted_prob", "truth", 
    "drugname_typaslab"), suffix = c(paste0(".", xlab), paste0(".", ylab)))
  
  p <- ggplot(joined, aes_string(x = paste0("prob.med.", xlab), y = paste0("prob.med.", ylab))) + 
    geom_point() + 
    facet_wrap( ~ predicted_prob, ncol = 2) + 
    labs(x = paste0("Prediction probs ", xlab), y = paste0("Prediction probs ", ylab), 
      title = title) + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  
  corrs <- by(joined, joined$predicted_prob, function(data) {
    cor(data[[paste0("prob.med.", xlab)]], data[[paste0("prob.med.", ylab)]], method = "spearman")
  })
  
  return(list(joined = joined, plot = p, spearman_corrs = corrs))
}

compare_onevsrest_models <- function(model1, model2, file) {
  # INPUT: 2 matrix_container_ext rows
  # OUTPUT: a pdf file comparing ROC and precision-recall curves and probabilities
  stopifnot(nrow(model1) == 1 & nrow(model2) == 1)
  # ROC curves:
  model1_ROC <- plot_perf_from_container(model1)
  model2_ROC <- plot_perf_from_container(model2)
  # precision-recall curves:
  model1_prc <- plot_perf_from_container(model1, what = "prec-recall")
  model2_prc <- plot_perf_from_container(model2, what = "prec-recall")
  # probabilities:
  # model1_probs <- plot_highest_probs(model1$pred_data[[1]], dosg_to_plot = "all", 
  #   which_moa = "protein_synthesis")
  # model2_probs <- plot_highest_probs(model2$pred_data[[1]], dosg_to_plot = "all", 
  #   which_moa = "protein_synthesis")
  
  model1_melted <- melt_pred_data(model1$pred_data[[1]])
  model2_melted <- melt_pred_data(model2$pred_data[[1]])
  prob_comparison <- compare_probabilities(model1_melted, model2_melted)$plot
  
  tmp_files <- c("A_tmp.pdf", "B_tmp.pdf")
  outdir <- c("./plots/")
  
  # ROC curves and precision-recall curves:
  pdf(file.path(outdir, tmp_files[1]), width = 12, height = 14)
  print(gridExtra::grid.arrange(model1_ROC, model2_ROC, model1_prc, model2_prc, nrow = 2))
  dev.off()
  
  # # prediction probabilities:
  # pdf(file.path(outdir, tmp_files[2]), width = 16, height = 28)
  # print(gridExtra::grid.arrange(model1_probs, model2_probs, nrow = 1))
  # dev.off()
  
  # comparison of probabilities
  pdf(file.path(outdir, tmp_files[2]), width = 8, height = 8)
  print(prob_comparison)
  dev.off()
  
  sys_cmd <- (paste0("pdftk ", paste0(outdir, tmp_files, collapse = " "), " cat output ", file))
  stopifnot(length(sys_cmd) == 1)
  system(sys_cmd)
  
  unlink(file.path(outdir, tmp_files))
  message("Output successfully saved to", file)
}

compare_multiclass_models <- function(model1_mc, model2_mc, file) {
  # input: 2 multi-class models
  # output: comparison of confusion matrices
  models <- list(model1 = model1_mc, model2 = model2_mc)
  conf_mats <- purrr::map(models, ~ plot_wide_confmat(get_wide_confmat(.x)))
  calib_plots <- purrr::map(models, plot_prob_calib)
  melted_models <- purrr::map(models, melt_pred_data, model_type = "multiclass")
  ## TO DO: FIX: this doesn't work:
  # compare_probabilities(melted_models$model1, melted_models$model2)
  model1_melt <- melted_models$model1
  model2_melt <- melted_models$model2
  prob_compar <- compare_probabilities(model1_melt, model2_melt)
  
  # too lazy to implement the same routine as in compare_onevsrest_models, should abstract 
  # the steps there first
  pdf(file, width = 12, height = 7)
  print(gridExtra::grid.arrange(conf_mats[[1]], conf_mats[[2]], ncol = 2))
  print(gridExtra::grid.arrange(calib_plots[[1]], calib_plots[[2]], ncol = 2))
  print(prob_compar$plot)
  dev.off()
}


