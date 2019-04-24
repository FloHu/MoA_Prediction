get_perf_measures_mcl <- function(resobj, perf_measures = list(mmce, kappa)) {
  # extract performance measures from a matrix container row
  # read: for each repetition of the nested CV: for each split of the repetition: 
  # get the performance measures
  perf_dfr <- imap_dfr(resobj, function(.outer_data, .outer_name) {
    imap_dfr(.outer_data, function(.inner_data, .inner_name) {
      perfs <- performance(.inner_data[["prediction"]], measures = perf_measures)
      df <- data.frame(cvrep = .outer_name, split = .inner_name, 
        stringsAsFactors = FALSE)
      df[, names(perfs)] <- perfs
      return(df)
    })
  })
  return(as_tibble(perf_dfr))
}

get_pred_data_mcl <- function(resobj) {
  # extract more detailed information on predictions from a matrix container row
  # input: results object
  pred_dfr <- imap_dfr(resobj, function(.outer_data, .outer_name) {
    imap_dfr(.outer_data, function(.inner_data, .inner_name) {
      pred_data <- .inner_data[["prediction"]][["data"]]
      pred_data[, c("cvrep", "split")] = list(.outer_name, .inner_name)
      pred_data <- as_tibble(pred_data)
      return(pred_data)
    })
  })
  # prettify output
  pred_dfr <- select(pred_dfr, drugname_typaslab, conc, cvrep, split, 
    prob.cell_wall:prob.protein_synthesis, response, truth) %>%
    arrange(cvrep, split, drugname_typaslab, conc)
  
  return(pred_dfr)
}

get_opt_path <- function(resobj) {
  # get the optimal hyperparameter path
  # input: results object
  # output: a data frame recording the performance for each hyperparameter 
  # combination (partial dependence was requested)
  opt_path_dfr <- 
    imap_dfr(resobj, function(.outer_data, .outer_name) {
      imap_dfr(.outer_data, function(.inner_data, .inner_name) {
        opt_path <- .inner_data[["opt_path"]][["data"]]
        opt_path[, c("cvrep", "split")] <- list(.outer_name, .inner_name)
        return(as_tibble(opt_path))
      })
    })
  
  return(opt_path_dfr)
}

get_opt_pars <- function(resobj) {
  # get the optimal hyperparameters 
  # input: results object
  # output: a data frame with the optimal hyperparameter set selected and the 
  # corresponding cross-validated performance
  dfr <- 
    imap_dfr(resobj, function(.outer_data, .outer_name) {
      imap_dfr(.outer_data, function(.inner_data, .inner_name) {
        opt_pars <- .inner_data[["opt_pars"]]
        opt_pars_perf <- .inner_data[["opt_pars_perf"]]
        opt_pars_df <- tibble(opt_pars = list(opt_pars), 
          mmce.cv.inner = opt_pars_perf["mmce.test.mean"], 
          kappa.cv.inner = opt_pars_perf["kappa.test.mean"])
        opt_pars_df[, c("cvrep", "split")] <- list(.outer_name, .inner_name)
        return(as_tibble(opt_pars_df))
      })
    })
  
  return(dfr)
}

get_pred_stabs <- function(pred_data) {
  # take prediction data ("pred_data") and get for each drug the sd of the 
  # probabilities for which it is predicted to be a particular moa
  # calculate sd for each drug-concentration combination separately 
  # since each moa gets one sd we will actually calculate the average sd 
  # (so average of 4 standard deviations)
  avg_sd_by_moa <- 
    group_by(pred_data, drugname_typaslab, conc) %>%
    summarise(sd_of_probs = mean(c(sd(prob.cell_wall), sd(prob.dna), 
      sd(prob.membrane_stress), sd(prob.protein_synthesis))))
  return(avg_sd_by_moa)
}

## alternative way of designing the above functions:
# 1) have a function operator make_resobj_accessor that modifies an 
# extractorfunction so that we don't have to keep recoding the mapping routine 
# problem: has to depend on a global "resobj" if we don't want the resulting 
# closures to become too large

# make_resobj_accessor <- function(extractorfunc) {
#   force(extractorfunc)
#   function(...) {
#     dfr <- 
#       imap_dfr(resobj, function(.x, .i) {
#         imap_dfr(.x, function(.xx, .ii) {
#           extractorfunc(.outer_data = .x, .outer_name = .i, .inner_data = .xx, 
#             .inner_name = .ii, ...)
#         })
#       })
#     return(dfr)
#   }
# }
# 
# perf_extractor_mcl <- function(.outer_data, 
#                                .outer_name, 
#                                .inner_data, 
#                                .inner_name, 
#                                perf_measures = list(mmce, kappa)) {
#   perfs <- performance(.inner_data[["prediction"]], measures = perf_measures)
#   df <- data.frame(cvrep = .outer_name, split = .inner_name, 
#     stringsAsFactors = FALSE)
#   df[, names(perfs)] <- perfs
#   return(as_tibble(df))
# }

# compare with above
# resobj <- res_new
# g <- make_resobj_accessor(extractorfunc = perf_extractor_mcl)
# g(perf_measures = list(mmce, kappa))
# output equivalent to: get_perf_measures_mcl(resobj = res_new)
# pryr::object_size(g)

# did not follow this route because sometimes we want to still do stuff with the 
# complete data frame and then stuff becomes difficult
