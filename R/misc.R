get_innerCV = function(repeat_ind, data = Rep_Nest_CV_instance){
  return(data[[repeat_ind]]$inner)
}

get_outerCV = function(repeat_ind, data = Rep_Nest_CV_instance){
  return(data[[repeat_ind]]$outer)
}

myhead <- function (x, n = 6L, ncol = 6L, rhs = FALSE) 
{
  # simply a modified version of head.data.frame showing just the first
  # ncol columns
  # argument rhs: if one wants to see the six columns from the right side
  stopifnot(length(n) == 1L)
  stopifnot(length(ncol) == 1L & (ncol > 0))
  n <- if (n < 0L)
    max(nrow(x) + n, 0L)
  else min(n, nrow(x))
  if (!rhs) {
    x[seq_len(n), seq_len(ncol), drop = FALSE]
  } else {
    from <- max(c(ncol(x) - ncol + 1, 1))
    x[seq_len(n), from:(ncol(x))]
  }
}

plot_pairplots <- function(genes, data = all_by_all) {
  # a very unflexible function that takes the all_by_all data frame and a vector with gene names 
  # to produce all possible combinations
  genes <- enquo(genes)
  
  filter(data, gene_synonym.x %in% !!genes, gene_synonym.y %in% !!genes) %>%
    ggplot(aes(x = s_score.x, y = s_score.y)) + 
    geom_point(aes(colour = process_broad), alpha = 0.5) + 
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    facet_wrap(gene_synonym.x ~ gene_synonym.y, strip.position = "bottom", scales = "free") + 
    theme_bw() + 
    labs(title = "Top strip = x-axis, bottom strip = y-axis")
}

filter_container <- function(container_obj, 
  dosgs = unique(container_obj$drug_dosages), 
  feats = unique(container_obj$feat_preselect), 
  chemfeats = NULL, 
  models = unique(container_obj$fitted_model)) {
  # this doesn't work:
  # print(as.list(match.call(expand.dots = FALSE)))
  # use this instead:
  # mget(names(formals()), sys.frame(sys.nframe()))
  my_args <- mget(names(formals())[-1], rlang::current_env())
  container_obj <- 
    container_obj %>%
    dplyr::filter(drug_dosages %in% my_args$dosgs, 
      feat_preselect %in% my_args$feats, 
      fitted_model %in% my_args$models, 
      if (is.null(my_args$chemfeats)) {
        TRUE
      } else {
        chemical_feats == my_args$chemfeats
      })
  
  return(container_obj)
}

get_outlier_boundaries <- function(data_vec) {
  # takes a numeric vector and returns the outlier boundaries according to boxplot.stats
  stopifnot(is.numeric(data_vec))
  outliers <- sort(boxplot.stats(data_vec)$out)
  lower_bound <- NA
  upper_bound <- NA
  if (any(outliers < 0)) {
    lower_bound <- max(outliers[outliers < 0])
  }
  if (any(outliers > 0)) {
    upper_bound <- min(outliers[outliers > 0])
  }
  
  data_vec_noout <- 
    if (is.na(lower_bound) & is.na(upper_bound)) { # no outliers = no boundaries
      data_vec
    } else if (!(is.na(lower_bound) | is.na(upper_bound))) { # i.e. we have 2 boundaries
      data_vec[data_vec > lower_bound & data_vec < upper_bound]
    } else if (is.na(lower_bound)) {
      data_vec[data_vec < upper_bound] 
    } else {
      data_vec[data_vec > lower_bound]
    }
  
  return(list(lower_bound = lower_bound, 
    upper_bound = upper_bound, 
    input_data = data_vec, 
    input_data_noout = data_vec_noout))
}

make_filename <- function(matrix_container_row) {
  filename <- paste0(paste(unlist(matrix_container_row[, c("hyperparam_grid_name", "drug_dosages", 
    "feat_preselect", "chemical_feats")]), collapse = "_"), ".rds")
  return(filename)
}

get_resobj_from_row <- function(matrix_container_row, dir) {
  # takes a matrix container row and reads in the results object
  # depends on function make_filename()
  return(readRDS(file.path(dir, make_filename(matrix_container_row))))
}

make_blocks <- function(dfr, blockvar) {
  # takes a data frame where observations belonging together are next to each 
  # other in the column blockvar, then returns a factor containing the blocks 
  # (for mlr)
  rl <- rle(dfr[[blockvar]])
  blocks <- factor(rep(seq_along(rl$lengths), rl$lengths))
  return(blocks)
}

nunique <- function(...) {
  length(unique(...))
}

sort_unique <- function(...) {
  sort(unique(...))
}

