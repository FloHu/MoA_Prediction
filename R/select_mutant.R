select_mutant <- function(dfr) {
  if (is.null(dfr$strain)) {
    stop("Did not find column 'strain'")
  }
  if ( length(unique(dfr$strain != 1)) ) {
    dfr <- group_by(dfr, strain) %>%
      mutate(n_signif_strain = sum(significant)) %>%
      ungroup() %>%
      filter(n_signif_strain == max(n_signif_strain)) %>%
      select(-n_signif_strain)
    # it's still possible that two strains happen to have same number of significant interactions
    # then we pick the first strain by alphabetical order
    if ( length(unique(dfr$strain != 1)) ) {
      first <- sort(dfr$strain)[1]
      dfr <- dfr[dfr$strain == first, ]
    }
    return(dfr)
  } else {
    return(dfr)
  }
}
