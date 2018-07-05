

corr_feat_selection = function(data, thres = NULL){
    corr_mat = cor(data)
    
    ut <- upper.tri(corr_mat)
    
    pair_corr = data.frame(featA = rownames(corr_mat)[row(corr_mat)[ut]],
                        featB = rownames(corr_mat)[col(corr_mat)[ut]],
                        cor  = (corr_mat)[ut], stringsAsFactors = FALSE)
    
    
    if(is.null(thres)){
        corr_values = sort(abs(pair_corr$cor))
        thres = findElbow(y = corr_values)
        thres = corr_values[thres]
    }
    
    stop = FALSE
    while(!stop){
        
        max = which.max(pair_corr$cor)
        # WARNING : threshold find with abs values, maybe should also work with this here
        if(pair_corr[max, "cor"] <= thres){
            stop = TRUE
        }
        
        # Put this output in a file
        print(pair_corr[max, ])
        # Put this output in a file
        
        
        a = pair_corr[max, "featA"]
        b = pair_corr[max, "featB"]
        a_corr = pair_corr %>% filter(featA == a | featB == a) %>% select(cor) %>% t() %>% mean() 
        b_corr = pair_corr %>% filter(featA == b | featB == b) %>% select(cor) %>% t() %>% mean() 
        
        if(a_corr >= b_corr){
            pair_corr = pair_corr %>% filter(featA != a & featB != a)
        }else{
            pair_corr = pair_corr %>% filter(featA != b & featB != b)
        }
    }
    
    toKeep = unique(c(pair_corr$featA, pair_corr$featB))
    return(data[, toKeep])
}

