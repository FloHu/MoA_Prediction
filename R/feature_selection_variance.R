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
