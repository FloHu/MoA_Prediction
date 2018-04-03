select_from_pca <- function(pca_obj, ncomps, feats_per_comp) {
   # takes from each of the first ncomps principal components the feats_per_comp features
   # that contribute most to the principal component loadings vector
   # if one of the features was already chosen in one of the earlier PCs, takes the next
   # feature until feats_per_comp unique features are identified per PC
   # pca_obj = principal component analysis object
   # ncomps = number of components that should be considered, starting with PC1
   # feats_per_comp: how many features should be picked per principal component
   if (ncomps < 1 |
       feats_per_comp < 1 |
       (! is.wholenumber(ncomps)) |
       (! is.wholenumber(feats_per_comp))) {
      stop("Must provide positive integer for ncomps and feats_per_comp")
   }
   n_features <- ncomps * feats_per_comp
   features <- vector(mode = "character", length = n_features) # where we collect all features
   m <- pca_obj$rotation
   features <- names(sort(abs(m[, 1]), decreasing = T)[1:feats_per_comp]) # ~ initialisation
   if (ncomps > 1) {
      for (pc in seq(2, ncomps)) {
         pool <- names(sort(abs(m[, pc]), decreasing = T))
         from <- 1
         to <- feats_per_comp
         candidates <- pool[c(from:to)]
         overlap <- intersect(features, candidates)
         while(length(overlap) > 0) {
            # if there is overlap, we need to advance in the feature space until no overlap is found
            candidates <- candidates[!candidates %in% overlap]
            from <- to + 1
            to <- to + length(overlap)
            if (to > length(pool)) {
               stop("Number of features exhausted, aborting. Try reducing feats_per_comp.")
            }
            candidates <- c(candidates, pool[c(from:to)])
            overlap <- intersect(features, candidates)
         }
         features <- c(features, candidates)
      }
   }
   return(features)
}
