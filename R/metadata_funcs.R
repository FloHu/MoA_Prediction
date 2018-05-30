# functions to add metadata to a matrix so that we can distinguish between:
# the feature_selection and transformation_procedure (if any) [top5pct, top10pct, ...]
# the dosage_selection procedure [most_interactions, ...]
# the datasets_included: currently one of nichols_2011, knime; to be extended later on
annotate_list <-
  function(l, add_metadata_func, ...) {
    # take a list of objects to which to apply one of the metadata functions below
    # ... takes the additional arguments for add_metadata_func
    l <- lapply(l, add_metadata_func, ...)
    return(l)
  }

add_metadata_featureselection <-
  function(drug_feature_matrix, feature_selection = c("all", "top5pct", "top10pct", "top15pct", "top20pct", "top25pct", "top30pct", "top40pct", "top50pct"))
  {
    drug_feature_matrix <- structure(drug_feature_matrix,
                                     feature_selection = match.arg(feature_selection))
    return(drug_feature_matrix)
  }

add_metadata_dosage <-
  function(drug_feature_matrix, dosage_selection = c("none", "most_interactions"))
  {
    drug_feature_matrix <- structure(drug_feature_matrix,
                                     dosage_selection = match.arg(dosage_selection))
    return(drug_feature_matrix)
  }

add_metadata_included_datasets <-
function(drug_feature_matrix, datasets_included = c("nichols_2011", "knime"))
  {
    drug_feature_matrix <- structure(drug_feature_matrix,
                                     datasets_included = match.arg(datasets_included, several.ok = TRUE))
    return(drug_feature_matrix)
  }

show_metadata <- function(drug_feature_matrix) {
  # takes a drug_feature_matrix and shows its attributes except for 'class', 'row.names' and 'names'
   attrs <- attributes(drug_feature_matrix)
   attrs[setdiff(names(attrs), c("class", "row.names", "names", "dim"))]
}
