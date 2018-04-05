features_selection_anova <- function(data , featStart = 3, featStop = ncol(data), pValueThreshold = 0.05, groupColName = "process_broad"){
    #Filter features in a data frame based on their capacity to seperate at least one group/factor level (ANOVA test)
    #Inputs
    #Data : raw data frame 
    #featStart : allow function to skip some column that are no predictive features, like the outcome variable
    toKeep <- apply(data[ ,featStart:featStop ], 2, function(x){
        a = aov(x ~ as.factor(data[ ,groupColName]))
        p = summary(a)[[1]][["Pr(>F)"]][[1]] #with Kruskal Wallis it is way easier to get the pValue
        (p <= pValueThreshold) #Return value
    }
    )
    
    toKeep = c(rep(TRUE, featStart-1), toKeep)
    return(data[ , toKeep])
    #possible improvement, add a shapiro test and perform AOV or KW based on the distribution assumption
    #But different tests, multiple pValue thresholds and pvalue correction won't have the same impact in both cases
}
