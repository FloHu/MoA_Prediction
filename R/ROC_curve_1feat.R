





ROC_curve_1feat = function(dt_mat, feat, moa){

    match.arg(arg = moa, choices = c("cell_wall", "dna", "membrane_stress", "protein_synthesis"))
    
    for(col_suppr in c("drugname_typaslab", "conc")){
        if(col_suppr %in% colnames(dt_mat)){
            dt_mat = dt_mat %>% select(-col_suppr)
        }
    }
    
    dt_mat = dt_mat[ ,c(feat, "process_broad") ] 
    threshold = sort(dt_mat[[feat]])

    rates <- lapply(threshold, FUN = function(t){
        
        pred_moa = dt_mat[which(dt_mat[, feat] <  t), ] %>% select(process_broad) %>% table()
        pred_not_moa = dt_mat[which(dt_mat[, feat] > t), ] %>% select(process_broad) %>% table()
            
        TP = ifelse(is.na(pred_moa[moa]), 0, pred_moa[moa])
        FN = ifelse(is.na(pred_not_moa[moa]), 0, pred_not_moa[moa])
            
        FP = sum(pred_moa) - TP
        TN = sum(pred_not_moa) - FN

        return(data.frame(tpr = TP/(TP+FN), fpr = FP/(FP+TN), row.names = NULL))
    })
    
    roc_data = data.frame(fpr = unlist(lapply(rates, function(x){x$fpr})), 
                          tpr = unlist(lapply(rates, function(x){x$tpr}))
                        )
   
    area = round(abs(computeAUC(roc_data)), digits = 3)
    
    p <- ggplot(data = roc_data, aes(x = fpr, y = tpr)) + geom_line(size = 1.5 ) +
        geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5) +
        labs(x = "False positive rate", y = "True positive rate",  title = paste(" ROC curve for feat :", feat,"impact on MoA :", moa, " AUC =", area )) +
        theme_bw()
    
    return(list(plot = p, data = roc_data, auc = area))
}




