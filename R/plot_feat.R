
plot_feat = function(RNCV_result = result_BT_10pc, ...){
    n_rep = length(RNCV_result)
    n_folds = length(RNCV_result[[1]])
    modelTree = T
    if(length(RNCV_result[[1]][[1]]) == 3){
        modelTree = F
    }
    
    all_rep_res = list()
    
    for (i in 1:n_rep){
        rep_res = list()
        for (j in 1:n_folds) {
            if(modelTree){
                rep_res[[j]] = getFeatureImportance(RNCV_result[[i]][[j]]$model)$res
            }else{
                rep_res[[j]] = RNCV_result[[i]][[j]][["featuresImportance"]]$res   
            }
        }
        all_rep_res[[i]] = rep_res
    }
    
    a = unlist(all_rep_res, recursive = F)
    aa = do.call(rbind.data.frame, a)
    aa = t(aa)

    
    colPal = colorRampPalette(c("white", "yellow","red"),space="rgb") 
    par(oma = c(0,2, 2,1))
    heatmap(aa, scale="column", col = colPal(10), ...)
    
    #Find cutoff
    #hist(sort(apply(aa, 1, sum), decreasing = T))
    #ab = aa[apply(aa, 1, sum) > 40, ]
    #heatmap(ab,  scale="none", col = colPal(10))    
    
    #library(plotly)
    #plot_ly(z = aa, colors = colPal(10), y = rownames(aa), type = "heatmap")
    
}
