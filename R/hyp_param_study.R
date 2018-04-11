
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
    return(a)   
}
