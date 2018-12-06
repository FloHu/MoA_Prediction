# by FH
get_resobj_from_row <- function(matrix_container_row, dir) {
  # takes a matrix container row and reads in the results object
  # depends on function make_filename()
  return(readRDS(file.path(dir, make_filename(matrix_container_row))))
}

get_feat_imps <- function(res_obj) {
  UseMethod("get_feat_imps")
}

get_feat_imps.ResampleResult <- function(res_obj) {
  feat_imps <- purrr::imap_dfr(res_obj$models, function(val, pos) {
    f <- getFeatureImportance(val)$res
    f <- tibble(gene = names(f), importance = unlist(f), modelnum = pos)
  })
  class(feat_imps) <- c("feat_imps_multiclass", class(feat_imps))
  return(feat_imps)
}

get_feat_imps.list <- function(res_obj) {
  # takes a result object and returns a tibble of feature importances
  all_combs <-
    list(cvrep = names(res_obj), fold = names(res_obj$`Nested CV 1`), moa = moas)
  all_combs <- as_tibble(expand.grid(all_combs, stringsAsFactors = FALSE))

  feat_imps <- pmap_dfr(all_combs, function(cvrep, fold, moa) {
    f <- getFeatureImportance(res_obj[[cvrep]][[fold]][[paste0("model_", moa)]])$res
    f <- tibble(gene = names(f), importance = unlist(f), cvrep = cvrep, fold = fold, moa = moa)
    return(f)
  })
  class(feat_imps) <- c("feat_imps_onevsrest", class(feat_imps))
  return(feat_imps)
}

summarise_feat_imps <- function(feat_imps) {
  UseMethod("summarise_feat_imps")
}

summarise_feat_imps.feat_imps_onevsrest <- function(feat_imps) {
  # get_feat_imps() returns a tibble where feature importance is recorded for
  # each feature for each cvrep + fold
  # this function summarises this information and adds columns
  # median_importance and 'pct_models_bin', indicating in how many % of models
  # a feature was used
  max_nmodels <- length(unique(feat_imps$cvrep)) * length(unique(feat_imps$fold))
  feat_imps_nest <-
    group_by(feat_imps, moa, gene) %>%
    # remove rows where feature importance = 0
    filter(importance != 0) %>%
    nest(.key = "importances") %>%
    mutate(pct_models = map_dbl(importances, ~ length(.x[["importance"]]) / max_nmodels),
      n_models = map_dbl(importances, ~ length(.x[["importance"]])),
      median_importance = map_dbl(importances, ~ median(.x[["importance"]]))) %>%
    select(-one_of("importances")) %>%
    arrange(moa, desc(pct_models), desc(median_importance))

    feat_imps_nest$pct_models_bin <- cut(feat_imps_nest$pct_models,
    breaks = c(0, 0.75, 0.8, 0.85, 0.9, 0.95, 1), labels = c("0-75%", "75-80%", "80-85%",
      "85-90%", "90-95%", "95-100%"))

    return(feat_imps_nest)
}

summarise_feat_imps.feat_imps_multiclass <- function(feat_imps) {
  max_nmodels <- length(unique(feat_imps$modelnum))
  feat_imps <- 
    group_by(feat_imps, gene) %>%
    filter(importance != 0) %>%
    nest(.key = "importances") %>%
    mutate(pct_models = map_dbl(importances, ~length(.x[["importance"]]) / max_nmodels), 
      n_models = map_dbl(importances, ~ length(.x[["importance"]])), 
      median_importance = map_dbl(importances, ~ median(.x[["importance"]]))) %>%
    select(-one_of("importances")) %>%
    arrange(desc(pct_models), desc(median_importance))
  
  feat_imps$pct_models_bin <- cut(feat_imps$pct_models,
    breaks = c(0, 0.75, 0.8, 0.85, 0.9, 0.95, 1), labels = c("0-75%", "75-80%", "80-85%",
      "85-90%", "90-95%", "95-100%"))
  return(feat_imps)
}

get_fingerprint <- function(feat_imps_s, log2thresh = NULL, top_n = NULL, method = c("thresh", "top_n")) {
  # takes a data frame with median feature importance across all models
  # user can choose to take all features above a certain threshold (after the importances have 
  # been log2 transformed) or the top_n features per MoA
  method <- match.arg(method)
  if (method == "thresh") {
    if (is.null(log2thresh)) {stop("Need to provide a threshold if using method 'thresh'")}
    feat_imps_s$log2median_importance <- log2(feat_imps_s$median_importance)
    fingerprint <- 
      dplyr::filter(feat_imps_s, log2median_importance > log2thresh) %>% 
      dplyr::arrange(gene)
  } else {
    if (is.null(top_n)) {stop("Need to provide top_n featurs per MoA if using method 'top_n'")}
    fingerprint <- 
      dplyr::group_by(fingerprint, moa) %>%
      dplyr::arrange(median_importance, .by_group = TRUE) %>%
      top_n(top_n, desc(median_importance))
  }
  return(fingerprint)
}

