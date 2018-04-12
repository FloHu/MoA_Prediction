
hyp_param_study = function(run_res){
    
    n_rep = length(run_res)
    n_folds = length(run_res[[1]])
    
    all_rep_res = list()
    
    for (i in 1:n_rep){
        rep_res = list()
        for (j in 1:n_folds) {
            rep_res[[j]] = run_res[[i]][[j]]$model$learner$par.vals
        }
        all_rep_res[[i]] = rep_res
    }

    a = unlist(all_rep_res, recursive = F)
    a = do.call(rbind.data.frame, a)
    
    #to test
    #a = cbind(a, sample(1:5, nrow(a), replace = T) )
    
    
    hypp_unique = apply(a, 2, function(x) length(unique(x))) > 1
    x = table(a[ , hypp_unique]) / length(a[ , hypp_unique])*100
    if(sum(hypp_unique) == 1){
        
        barplot(x, main = names(hypp_unique[hypp_unique]), col = rainbow(length(x)))
    }else(
        
        mosaicplot(x, color = rainbow(ncol(x)))
    )
    
    
    
    return(a)   
}
