corr_feat_selection = function(data, thres = NULL){
    
    if(!require(Hmisc)){
        install.packages("Hmisc")
    }
    library(Hmisc)


    corr_mat = rcorr(as.matrix(data), type="pearson")
    ut = upper.tri(corr_mat$r)
    
    pair_corr = data.frame(featA = rownames(corr_mat$r)[row(corr_mat$r)[ut]],
                        featB = rownames(corr_mat$r)[col(corr_mat$r)[ut]],
                        cor  = (corr_mat$r)[ut],
                        pval  = (corr_mat$P)[ut], stringsAsFactors = FALSE)
    #pair_corr= as.tbl(pair_corr)
    pair_corr$cor = abs(pair_corr$cor)
   

    if(is.null(thres)){
        corr_values = sort(abs(pair_corr$cor))
        thres = findElbow(y = corr_values)
        thres = corr_values[thres]
    }else if(thres == "bonf01"){
        pair_corr$bonf = p.adjust(pair_corr$pval, method ="bonferroni")
        thres = min(abs(pair_corr %>% filter(bonf <= 0.01) %>% select(cor)))
    }else if(thres == "bonf05"){
        pair_corr$bonf = p.adjust(pair_corr$pval, method ="bonferroni")
        thres = min(abs(pair_corr %>% filter(bonf <= 0.05) %>% select(cor)))
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
        
        a = pair_corr[max, ]$featA
        b = pair_corr[max, ]$featB
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

# To get a matrix of corr coef AND p values :
# ipak(Hmisc)
#pair_corr <- rcorr(as.matrix(data), type="pearson")
