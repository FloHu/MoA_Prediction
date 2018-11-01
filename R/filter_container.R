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
