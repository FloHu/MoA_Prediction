retrieve_medoids <- function(clusts, distm) {
   # clusts = output of cutree
   # distm: matrix to get distances between cluster members from
   # returns: named vector where the names are the medoids of cluster
   clusterlabels <- unique(clusts)
   medoids <- vector(mode = "numeric", length = length(clusterlabels))
   names(medoids) <- vector(mode = "character", length = length(clusterlabels))
   for (c in clusterlabels){
      clustermembers <- clusts[clusts == c]
      if (length(clustermembers) > 1) {
         for (gene in names(clustermembers)) {
            clustermembers[gene] <- mean(distm[gene, names(clustermembers)[ !names(clustermembers) %in% gene]])
         }
         medoids[c] <- min(clustermembers)
         names(medoids)[c] <- names(clustermembers)[which.min(clustermembers)]
      } else {
         medoids[c] <- clustermembers[1]
         names(medoids)[c] <- names(clustermembers[1])
      }
   }
   return(names(medoids))
}
