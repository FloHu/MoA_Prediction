feature_selection_variance <- function(data , featStart = 3, featStop = ncol(data), percentTop = 20){
  #Filter features in a data frame based on their variance. Keep only the top n% (by default 20)
  #Inputs
  #Data : raw data frame
  #featStart : allow function to skip some column that are no predictive features, like the outcome variable
  allVar <- apply(data[ ,featStart:featStop], 2, var)
  allVar <- sort(allVar, decreasing = T)
  
  allVar = allVar[1:(length(allVar)*(percentTop/100))]
  toKeep = c(names(data)[1:(featStart - 1)], names(allVar))
  # sanity check
  stopifnot(unique(toKeep) == toKeep)
  return(data[ , toKeep])
}

select_dosage_most_ias <- function(dfr) {
  subfr <- select_if(ungroup(dfr), is.logical) # need to ungroup!
  sum_signif <- apply(subfr, 1, sum)
  # first occurrence wins, meaning the highest dosage will have the maximum rank in the case of ties:
  # (as dosages were sorted previously)
  ranks <- rank(sum_signif, ties.method = "first")
  ranks <- ranks == max(ranks)
  dfr$selector <- ranks
  return(dfr)
}

select_mutant <- function(dfr) {
  if (is.null(dfr$strain)) {
    stop("Did not find column 'strain'")
  }
  if (length(unique(dfr$strain)) != 1) {
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

