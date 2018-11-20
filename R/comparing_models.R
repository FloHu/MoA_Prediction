melt_pred_data <- function(pred_data, model_type = c("onevsrest", "multiclass")) {
  # input: prediction data frame, either from slot $pred$data from mlr prediction objects or 
  # from matrix_container_ext
  model_type <- match.arg(model_type)
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
  
  joined <- inner_join(melted_multi, melted_onevsrest, by = c("conc", "predicted_prob", "truth", 
    "drugname_typaslab"), suffix = c(paste0(".", xlab), paste0(".", ylab)))
  
  ggplot(joined, aes_string(x = paste0("prob.med.", xlab), y = paste0("prob.med.", ylab))) + 
    geom_point() + 
    facet_wrap( ~ truth, ncol = 2) + 
    labs(x = paste0("Prediction probs ", xlab), y = paste0("Prediction probs ", ylab), 
      title = title)
}

(melted_multi <- melt_pred_data(resampled_multiclass_pred, model_type = "multiclass"))
(melted_onevsrest <- melt_pred_data(matrix_container_ext$pred_data[[the_line]], 
  model_type = "onevsrest"))

compare_probabilities(a = melted_onevsrest, b = melted_multi, title = "Prob comparison")



