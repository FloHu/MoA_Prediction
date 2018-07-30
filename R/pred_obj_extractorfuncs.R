# extractor functions get the predictions/ThreshVsPerfData from each test fold of each prediction
# object of each repat of the cross-validation
# the actual extraction is performed in prediction_merger()
# Measures as an argument with default values ? Thus same function for ROC or prec racal or mcc, etc VS threshold curves

perf_extractor <- function(x, moa) {
   # x = a prediction object
   # moa = mode of action
   perfs <- generateThreshVsPerfData(x[[paste0("prediction_", moa)]], measures = list(fpr, tpr, ppv))$data
   return(perfs)
}

pred_extractor <- function(x, moa) {
   # x = a prediction object
   # moa = mode of action
   x <- as.data.frame(x[[paste0("prediction_", moa)]])
   toreplace <- c(paste0("prob.", moa), paste0("prob.not_", moa))
   indstoreplace <- match(toreplace, names(x))
   names(x)[indstoreplace] <- c("prob.moa", "prob.not_moa")
   return(x)
}

prediction_merger <- function(resultsobj, moa, extractorfunc) {
   # resultsobj = run from a repeated nested CV
   # moa = one of our moas - check
   # extractorfunc = function that pulls out the stuff from the prediction object
   dfr_list <-
   map(resultsobj, function(.x) { # apply to each repeat of the nested CV
      map2(.x, names(.x), function(.x, .y) { # apply to each fold, record the name
         dfr <- extractorfunc(.x, moa)
         dfr$fold <- .y
         return(dfr)
      }) %>%
      bind_rows() # to rbind all the test folds
   })

   # collapse into a data frame
   dfr <- imap_dfr(dfr_list, function(.x, .y) {
      .x[["cvrep"]] <- .y
      return(.x)
   })

   # add moa information
   dfr[["moa_modelled"]] <- moa
   return(dfr)
}
