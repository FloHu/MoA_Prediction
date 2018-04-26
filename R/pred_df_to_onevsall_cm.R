pred_df_to_onevsall_cm <- function(pred_df, targetclass, truthcol = "truth", predictcol = "response",
                                   probability_col = NULL, thresh = NULL) {
   # turns a data frame into a one (targetclass) vs all confusion matrix
   # can adjust the threshold for class prediction of targetclass by passing argument thresh
   # otherwise truth and response columns are not modified
   # must provide a column containing the probability values (probability_col argument)
   pred_df[, c(truthcol, predictcol)] <-
      lapply(pred_df[, c(truthcol, predictcol)], function(x) {
         x <- as.character(x)
         x <- ifelse(x == targetclass, targetclass, "other")
         return(x)
      })

   if (! is.null(thresh)) {
      if (is.null(probability_col)) stop("Must provide a column for probability thresholding")
      pred_df[[predictcol]] <- ifelse(pred_df[[probability_col]] > thresh, targetclass, "other")
   }

   pred_df[, c(truthcol, predictcol)] <-
      lapply(pred_df[, c(truthcol, predictcol)], function(x) {
         x <- factor(x, levels = c(targetclass, "other"))
         return(x)
      })

   cm <- table(pred_df[, c(truthcol)],
               pred_df[, c(predictcol)])
   names(dimnames(cm)) <- list(truthcol, predictcol)

   return(cm)
}
