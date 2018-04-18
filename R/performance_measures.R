calculate_prec <- function(m) {
   # expects a 2x2 matrix, will return precision for the first of the two classes
   if (! (all(dim(m) == c(2, 2)))) {
      stop("Function expects a 2x2 matrix")
   }

   precision <- m[1, 1] / colSums(m)[1]

   classname <- dimnames(m)[[1]][1]
   if (is.null(classname)) {
      classname <- "`1`"
   }
   names(precision) <- classname
   return(precision)
}

calculate_recall <- function(m) {
   # expects a 2x2 matrix, will return precision for the first of the two classes
   if (! (all(dim(m) == c(2, 2)))) {
      stop("Function expects a 2x2 matrix")
   }

   recall <- m[1, 1] / rowSums(m)[1]

   classname <- dimnames(m)[[1]][1]
   if (is.null(classname)) {
      classname <- "`1`"
   }
   names(recall) <- classname
   return(recall)
}

calculate_fpr <- function(m) {
   # expects a 2x2 matrix, will return precision for the first of the two classes
   if (! (all(dim(m) == c(2, 2)))) {
      stop("Function expects a 2x2 matrix")
   }

   fpr <- m[2, 1] / rowSums(m)[2]

   classname <- dimnames(m)[[1]][1]
   if (is.null(classname)) {
      classname <- "`1`"
   }
   names(fpr) <- classname
   return(fpr)
}
