


# Take prediction data from contaimer matrix
# plot drug probabilities repartition across repetitions
drugs_pred_prob_from_container = function(pred_data, moa, order = median, ...){

    pred_data = pred_data %>% filter(moa_modelled == moa) %>% select(prob.moa, drugname_typaslab, truth)
    mean_pred = by(data = pred_data$prob.moa, pred_data$drugname_typaslab, FUN = order)
    grp_pred = split(pred_data$prob.moa, f = pred_data$drugname_typaslab)
    a = names(sort(mean_pred))
    grp_pred = grp_pred[a]
    
    drug_moa = pred_data %>% filter(truth == moa) %>% select(drugname_typaslab) %>% unique()
    
    colorMap = rep("#FF0000FF", length(grp_pred))
    colorMap[which(a %in% drug_moa$drugname_typaslab)] = "#00FFFFFF"
    
    par(mar = c(5,10,4,2))
    boxplot(grp_pred, lwd = 1.5, col = colorMap,  las = 2, horizontal = TRUE, cex.axis = 0.8, ...)
    legend("bottomright", inset = 0.05, legend = c(moa, paste0("not_", moa )), fill = c("#00FFFFFF", "#FF0000FF"))
    return(a)
}


# Same but compare with/without chemFeat
prob_shift_plot = function(pred_data, pred_data_noChem, moa, order = median, ...){
    
    pred_data = pred_data %>% filter(moa_modelled == moa) %>% select(prob.moa, drugname_typaslab, truth)
    pred_data_noChem = pred_data_noChem %>% filter(moa_modelled == moa) %>% select(prob.moa, drugname_typaslab, truth)
    
    mean_pred = by(data = pred_data$prob.moa, pred_data$drugname_typaslab, FUN = median)
    grp_pred = split(pred_data$prob.moa, f = pred_data$drugname_typaslab)
    ord = names(sort(mean_pred))
    grp_pred = grp_pred[ord]
    
    par(mar = c(5,10,4,2))
    boxplot(grp_pred, lwd = 1,  las = 2, horizontal = TRUE, cex.axis = 0.8, col = "lightblue",  at = seq(from = 1, to = 3*length(grp_pred), by = 3), outline = F)
    
    grp_pred = split(pred_data_noChem$prob.moa, f = pred_data_noChem$drugname_typaslab)
    ord = names(sort(mean_pred))
    grp_pred = grp_pred[ord]
    
    boxplot(grp_pred, lwd = 1,  las = 2, horizontal = TRUE, cex.axis = 0.8, col = "red", add = T, at = seq(from = 2, to = 3*length(grp_pred), by = 3), yaxt = 'n',  outline = F, ...)
}



#Compare Dosages prediction
drugs_pred_prob_by_dosage = function(pred_data, moa, order = median){
    
    
    pred_data = pred_data %>% filter(moa_modelled == moa) %>% select(prob.moa, drugname_typaslab, truth, conc)
    # Redefine name in order to group by drugs and by dosage
    pred_data$conc = paste(pred_data$drugname_typaslab, pred_data$conc, sep = "_")
    
    mean_pred = by(data = pred_data$prob.moa, pred_data$conc, FUN = order)
    grp_pred = split(pred_data$prob.moa, f = pred_data$conc)
    ord = names(sort(mean_pred))
    grp_pred = grp_pred[ord]
    
    names(grp_pred) = lapply(str_split(names(grp_pred), pattern = "_"), function(x){x[1]})

    ggplot(data = pred_data, aes(x = conc, y = prob.moa )) + geom_boxplot() + 
        facet_wrap(~drugname_typaslab, scales = "free_y",nrow =18, ncol = 4 ) + 
        coord_flip() + 
        theme_bw()
    
    
}
