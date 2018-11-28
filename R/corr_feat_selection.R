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

