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
  models = unique(container_obj$fitted_model), 
  get_line_number = FALSE) {
  # this doesn't work:
  # print(as.list(match.call(expand.dots = FALSE)))
  # use this instead:
  # mget(names(formals()), sys.frame(sys.nframe()))
  my_args <- mget(names(formals())[-1], rlang::current_env())
  
  line_number <- with(container_obj, drug_dosages %in% my_args$dosgs & 
      feat_preselect %in% my_args$feats & 
      fitted_model %in% my_args$models & 
      if (is.null(my_args$chemfeats)) {TRUE} else {chemical_feats == my_args$chemfeats})
  
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
  
  if (get_line_number) 
    return(which(line_number))
  else 
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

fit_or_load <- function(varname, directory, dfm, fitting_fun, ...) {
  # function first checks if there is a file called varname in directory 
  # if not, runs the function fitting_fun on dfm (drug-feature matrix) 
  # and then saves the output to disc
  # ... = further arguments passed to fitting_fun
  my_filename <- file.path(directory, paste0(varname, ".rds"))
  if (file.exists(my_filename)) {
    message("Loaded ", paste0(varname, ".rds"), " from disk")
    return(readRDS(my_filename))
  } else {
    message("File not found, fitting model on ", deparse(substitute(dfm)), 
      " using ", deparse(substitute(fitting_fun)), "\n\n")
    my_fit <- fitting_fun(dfm, ...)
    saveRDS(my_fit, file = my_filename)
    message("Fitting successful, wrote file to ", my_filename)
    return(my_fit)
  }
}

generate_corpairs <- function() {
  # just outsourced for readability
  corpairs <- melt_cormat_to_pairs(feats_cor)
  
  # retrieving all protein complexes, code from notebook 1 (Leonard) -----------
  cmplx = read_tsv(file = "./data/All_instances_of_Protein-Complexes_in_Escherichia_coli_K-12_substr._MG1655.txt")
  (load("./data/genesWithEG_ID.RData"))
  gene_synonyms$synonym = toupper(gene_synonyms$synonym)
  
  genes_cmplx = lapply(cmplx$`Genes of polypeptide, complex, or RNA`,
    FUN = function(x){
      genenames = str_match_all(string = x, pattern = "[a-zA-Z]{3,4}")
      genenames = toupper(unlist(genenames))
    })
  
  names(genes_cmplx) = cmplx$`Protein-Complexes`
  genes_cmplx = genes_cmplx[!is.na(genes_cmplx)]
  
  lens = lapply(genes_cmplx, function(x){length(x)})
  genes_cmplx = genes_cmplx[lens != 1]
  saveRDS(genes_cmplx, file = file.path(outdir, "genes_cmplx.rds"))
  
  # add to each pair the information if the two genes are also part of the same 
  # protein complex
  genes_cmplx_v <- unique(unname(unlist(genes_cmplx)))
  
  corpairs$in_same_cmplx <- 
    map2_lgl(corpairs$featA, corpairs$featB, function(.featA, .featB) {
      drugpair <- c(.featA, .featB)
      if (!all(c(.featA, .featB) %in% genes_cmplx_v)) {
        FALSE
      } else {
        any(map_lgl(genes_cmplx, ~ (sum(.x %in% drugpair) == 2)))
      }
    })
  
  return(corpairs)
}

make_my_task <- function(dfm, targetvar = "process_broad", blockvar = NULL) {
  # dfm = drug-feature matrix
  # this function just removes drugname_typaslab, conc; then defines 
  # process_broad as target column for a task
  # can provide a blocking variable 
  dfm <- as.data.frame(dfm)
  if (!is.null(blockvar)) blocks <- make_blocks(dfm, blockvar = blockvar)
  
  to_remove <- c("drugname_typaslab", "conc")
  df <- dfm[, !colnames(dfm) %in% to_remove]
  
  le_task <- makeClassifTask(data = df, target = targetvar, 
    blocking = if (!is.null(blockvar)) blocks else NULL)
  
  # annoying thing about mlr task is that there is no option to save metadata 
  # about the task 
  # so we generate our own custom slot
  le_task$data_complete <- as_tibble(dfm)
  
  return(le_task)
}

melt_pred_data_mcl <- function(pred_data) {
  # input: pred_data from mc_ext
  # output: melted data frame that can be used to plot a prediction heatmap
  
  my_cols <- colnames(pred_data)
  
  if (!("conc" %in% my_cols)) {
    stop("Specified data frame doesn't contain a drug concentration")
  }
  
  # melt and average probabilities
  
  tryCatch(
    {melted <- tidyr::gather(pred_data, prob.cell_wall:prob.protein_synthesis, 
      key = "predicted_prob", value = "probability") %>%
      group_by(conc, predicted_prob, truth, drugname_typaslab) %>%
      summarise(prob.min = boxplot.stats(probability)$stats[2], 
        prob.max = boxplot.stats(probability)$stats[4], 
        prob.med = median(probability)) %>%
      ungroup()
    }, 
    error = function(cnd) {
      stop("It looks like the data frame doesn't contain the four main modes of action.")
    }
  )
  
  return(melted)
}

