melt_cormat_to_pairs <- function(cor_mat) {
  # input: a correlation matrix
  # output: a tibble with each row being one correlation pair
  if (!(is.matrix(cor_mat) && is.numeric(cor_mat) && (dim(cor_mat)[1L] == dim(cor_mat)[2L]))) {
    stop("Need to provide a numeric correlation matrix as input")
  }
  
  ut <- upper.tri(cor_mat)
  melted_cor_mat <- tibble::tibble(
    featA = rownames(cor_mat)[row(cor_mat)[ut]],
    featB = rownames(cor_mat)[col(cor_mat)[ut]],
    cor = cor_mat[ut]
  )
  
  return(melted_cor_mat)
}

merge_features <- function(dfr, cluster_fr, f) {
  # takes a data frame and a list of features
  # each element of the list indicates which columns of the data frame should 
  # be merged using function f
  stopifnot(all(unlist(cluster_fr$members) %in% colnames(dfr)))
  to_merge <- cluster_fr$members
  names(to_merge) <- cluster_fr$clust_label
  
  newfeats <- imap_dfc(to_merge, function(featureset, label) {
    v <- tibble(apply(dfr[, featureset], 1, f))
    names(v) <- label
    return(v)
  })
  
  message("Mapped ", length(unlist(to_merge)), " into ", nrow(cluster_fr), " derived features.")
  dfr <- select(dfr, -one_of(unlist(to_merge)))
  dfr <- bind_cols(dfr, newfeats)
  return(dfr)
}





