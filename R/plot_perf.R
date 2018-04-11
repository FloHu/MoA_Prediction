

plot_perf = function(RNCV_result = swag, measure = mmce, type = "all"){
    n_rep = length(RNCV_result)
    n_folds = length(RNCV_result[[1]])
    
    all_rep_res = list()
    
    for (i in 1:n_rep){
        rep_res = c()
        for (j in 1:n_folds) {
            rep_res = c(rep_res, performance(RNCV_result[[i]][[j]]$prediction, measures = list(measure)))
        }
        all_rep_res[[i]] = rep_res
    }
    
    
    if(type == "all"){
        boxplot(unlist(all_rep_res), lwd = 1.5)
        stripchart(unlist(all_rep_res), method = "jitter", pch = 21, col = "black", bg = "orange", cex = 2, vertical = T, add = T)
    }
    
    if(type == "byrep"){
        boxplot(all_rep_res)
        stripchart(all_rep_res, method = "jitter", pch = 21, col = "black", bg = "orange", cex = 2, vertical = T, add = T)
    }

}
