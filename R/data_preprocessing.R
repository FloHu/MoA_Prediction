feature_selection_variance <- function(data, 
                                       percentTop, 
                                       featStart = 3, 
                                       featStop = ncol(data)) {
  # Filter features in a data frame based on their variance. Keep only the top n%
  # Inputs:
  # data: raw data frame
  # featStart  allow function to skip some column that are not predictive 
  # features, like the outcome variable
  allVar <- apply(data[, featStart:featStop], 2, var)
  allVar <- sort(allVar, decreasing = T)
  
  allVar = allVar[1:(length(allVar) * (percentTop/100))]
  toKeep = c(names(data)[1:(featStart - 1)], names(allVar))
  # sanity check
  stopifnot(unique(toKeep) == toKeep)
  return(data[, toKeep])
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

select_most_ias <- function(the_matrix, the_matrix_signifs) {
  # purpose: select dosage with most interactions of the_matrix, if there are 
  # ties, take the highest dosage
  # INPUTS:
  #  * the_matrix: a drug-feature matrix with several concentrations per drug
  #  * the_matrix_signifs: a matrix containing at least the same drug-conc 
  #    combinations of the_matrix, the same features, but indicating whether 
  #    the interaction is significant or not
  rowsa <- paste(the_matrix$drugname_typaslab, the_matrix$conc)
  rowsb <- paste(the_matrix_signifs$drugname_typaslab, the_matrix_signifs$conc)
  colsa <- colnames(the_matrix)
  colsb <- colnames(the_matrix_signifs)
  metavars <- c("drugname_typaslab", "conc")
  
  stopifnot({
    all(colsa %in% colsb)
    all(metavars %in% colnames(the_matrix))
    all(rowsa %in% rowsb)
    all(sapply(select(the_matrix_signifs, -one_of(metavars)), is.logical))
  })
  
  the_matrix_signifs <- the_matrix_signifs[rowsb %in% rowsa, colsb %in% colsa]
  the_matrix_signifs <- bind_cols(the_matrix_signifs[, metavars], 
    nsignif = apply(select(the_matrix_signifs, -one_of(metavars)), 1, sum))
  
  winning_concs <- 
    left_join(the_matrix[, metavars], the_matrix_signifs) %>%
    group_by(drugname_typaslab) %>%
    arrange(drugname_typaslab, desc(nsignif), desc(conc)) %>%
    slice(1)
  
  the_matrix_subset <- semi_join(the_matrix, winning_concs) %>%
    arrange(drugname_typaslab)
  stopifnot({
    ncol(the_matrix) == ncol(the_matrix_subset)
    all(the_matrix$drugname_typaslab %in% the_matrix_subset$drugname_typaslab)
  })
  return(the_matrix_subset)
}

select_mutant <- function(dfr) {
  # DEPRECATED
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

