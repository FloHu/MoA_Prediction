add_metadata <- function(drug_feature_matrix,
                         drug_selection = c("not_specified", "all_drugs", "main_moa_drugs"),
                         drug_classification = c("not_specified", "before_2018-05-09", "after_2018-05-09"),
                         feature_selection = c("not_specified", "all", "top5pct", "top10pct", "top15pct", "top20pct", "top25pct", "top30pct", "top40pct", "top50pct"),
                         dosage_selection = c("not_specified", "most_interactions"),
                         datasets_included = NULL) {
   # adds metadata to a matrix so that we can distinguish between:
   # the drugs that were selected: attribute drug_selection [all, main_moa_drugs, ...]
   # and similarly for the drug_classification used [before_2018-05-09, after_2018-05-09, ...]
   # the feature_selection and transformation_procedure (if any) [top5pct, top10pct, ...]
   # the dosage_selection procedure [most_interactions, ...]
   # the datasets_included: currently one of nichols_2011, knime; to be extended later on
   # the meanings of these metadata are documented in the notebook
   drug_feature_matrix <- structure(drug_feature_matrix,
                                    drug_selection = match.arg(drug_selection),
                                    drug_classification = match.arg(drug_classification),
                                    feature_selection = match.arg(feature_selection),
                                    dosage_selection = match.arg(dosage_selection),
                                    datasets_included = match.arg(datasets_included, choices = c("not_specified", "nichols_2011", "knime"), several.ok = TRUE),
                                    date_created = Sys.Date())
   return(drug_feature_matrix)
}

show_metadata <- function(drug_feature_matrix) {
  # takes a drug_feature_matrix and shows its attributes except for 'class', 'row.names' and 'names'
   stopifnot(is.data.frame(drug_feature_matrix))
   attrs <- attributes(drug_feature_matrix)
   attrs[setdiff(names(attrs), c("class", "row.names", "names"))]
}
