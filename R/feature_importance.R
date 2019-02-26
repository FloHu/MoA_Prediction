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
  # med_imp and 'pct_models_bin', indicating in how many % of models
  # a feature was used
  max_nmodels <- length(unique(feat_imps$cvrep)) * length(unique(feat_imps$fold))
  feat_imps_nest <-
    group_by(feat_imps, moa, gene) %>%
    # remove rows where feature importance = 0
    filter(importance != 0) %>%
    nest(.key = "importances") %>%
    mutate(pct_models = map_dbl(importances, ~ length(.x[["importance"]]) / max_nmodels),
      n_models = map_dbl(importances, ~ length(.x[["importance"]])),
      med_imp = map_dbl(importances, ~ median(.x[["importance"]]))) %>%
    select(-one_of("importances")) %>%
    arrange(moa, desc(pct_models), desc(med_imp))

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
      med_imp = map_dbl(importances, ~ median(.x[["importance"]]))) %>%
    select(-one_of("importances")) %>%
    arrange(desc(pct_models), desc(med_imp))
  
  feat_imps$pct_models_bin <- cut(feat_imps$pct_models,
    breaks = c(0, 0.75, 0.8, 0.85, 0.9, 0.95, 1), labels = c("0-75%", "75-80%", "80-85%",
      "85-90%", "90-95%", "95-100%"))
  return(feat_imps)
}

get_fingerprint <- function(feat_imps_s, log2thresh = NULL, top_n = NULL) {
  # takes a data frame with median feature importance across all models
  # user can choose to take all features above a certain threshold (after the importances have 
  # been log2 transformed) - this doesn't consider the MoA 
  # if top_n is chosen, the top_n features per MoA are picked
  if (!xor(is.null(log2thresh), is.null(top_n))) {
    stop("Have to provide either log2thresh or top_n, not both")
  }
  if (is.null(top_n)) {
    feat_imps_s$log2med_imp <- log2(feat_imps_s$med_imp)
    fingerprint <- 
      dplyr::filter(feat_imps_s, log2med_imp > log2thresh) %>% 
      dplyr::arrange(gene)
  } else {
    fingerprint <- 
      dplyr::group_by(feat_imps_s, moa) %>%
      dplyr::arrange(desc(med_imp), .by_group = TRUE) %>%
      top_n(top_n, med_imp)
  }
  return(fingerprint)
}

plot_importance_distribs <- function(feat_imps, fingerprint, title) {
  UseMethod("plot_importance_distribs")
}

plot_importance_distribs.feat_imps_onevsrest <- function(feat_imps, fingerprint) {
    # takes a feat_imps_onevsrest tibble and plots distribution of feature 
    # importances across all nested CV repeats for each MoA
    # only the features in fingerprint are kept
    my_moas <- unique(feat_imps$moa)
    moa_cols <- c(cell_wall = "#1b9e77", dna = "#d95f02", 
      membrane_stress = "#7570b3", protein_synthesis = "#e7298a")
    
    to_plot <- 
      dplyr::filter(feat_imps, paste(gene, moa) %in% 
          paste(fingerprint$gene, fingerprint$moa)) %>%
      group_by(gene, moa, cvrep) %>%
      dplyr::summarise(med_imp = median(importance), 
        med_imp_log2 = log2(median(importance))) %>%
      ungroup()
    to_plot$colour <- moa_cols[to_plot$moa]
    
    imp_list <- list()
    imp_list[my_moas] <- list(NULL)
    imp_list <- imap(imp_list, function(.el, .name) {
      .el <- filter(to_plot, moa == .name) %>%
        mutate(gene = fct_reorder(gene, med_imp))
      p <- ggplot(.el, aes(x = gene, y = med_imp_log2, colour = colour)) + 
        stat_summary(geom = "pointrange", fun.y = "median", size = 0.5, 
          fun.ymin = function(x) {fivenum(x)[2]}, 
          fun.ymax = function(x) {fivenum(x)[4]}) + 
        facet_wrap( ~ moa, ncol = 2, scales = "free") + 
        coord_flip() + 
        labs(x = "Feature", y = "Feature importance (log2)") + 
        scale_colour_identity()
      return(list(imps = .el, plot = p))
    })
    
    return(imp_list)
  }

