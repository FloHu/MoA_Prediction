

#Concatenate prediction of test set from each fold of the outer CV
#Get prediction object for one model of 1vs ALL
#Each individual is predected once
cat_outer_fold_pred = function(res = result_xgb_5pc_4model, repetition = 1, moa = "dna"){
    
    if(!moa %in% c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
        print("Invalid Mode of action")
        return(-1)
    }
    
    n_fold = length(res[[1]])
    output = c()

    for(i in 1:n_fold){
        output = rbind(output, res[[repetition]][[i]][[paste0("prediction_", moa)]]$data)
    }

    predObj = res[[1]][[1]][[paste0("prediction_", moa)]]
    predObj$data = output
    
    return(predObj)
}

#Plot function work for 1 MoA or all
#Generate a boxplot of a measure across all repetition
#moa = "all" allows to compare all classes 
plot_mean_perf = function(res = result_xgb_5pc_4model, meas = auc, moa = "dna", ...){
    
    set_MoA = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    
    if(!moa %in% c(set_MoA, "all")){
        print("Invalid Mode of action")
        return(-1)
    }
    
    if(moa == "all"){
        toPlot = c()
        for (m in set_MoA){
            x = c()
            for(rep in 1:10){
                pred = cat_outer_fold_pred(res = res, repetition = rep, moa = m)
                x = c(x, performance(pred, measures = meas))
            }
            toPlot = rbind(toPlot, x)
            
        }
        rownames(toPlot) = set_MoA
        boxplot(t(toPlot), lwd = 1.5, outline = F, 
                main = paste0("Measure comparison across MoA :", as.character(meas$id) ),
                col = rainbow(length(set_MoA)),...)
    }else{
        toPlot = c()
        for(rep in 1:10){
            pred = cat_outer_fold_pred(repetition = rep, moa = moa)
            toPlot = c(toPlot, performance(pred, measures = meas))
        } 
        boxplot(toPlot, lwd = 1.5, outline = F, 
                main = paste0("Measure comparison across Repetition :", as.character(meas$id) ) ,
                sub = paste0("MoA : ", moa), ...)
        stripchart(toPlot, method = "jitter", pch = 21, col = "black", bg = rainbow(10), cex = 2, vertical = T, add = T, ...)
    }

}



#Plot function work for 1 MoA or all
#Generate a ROC curve of all repetition, each individual is then predicted 10 times
#moa = "all" allows to compare all classes 
plot_ROC_rep = function(res = result_xgb_5pc_4model, moa = "dna"){

    set_MoA = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    
    if(!moa %in% c(set_MoA, "all")){
        print("Invalid Mode of action")
        return(-1)
    }
    if(moa != "all"){
        all_rep = cat_outer_fold_pred(res = res, repetition = 1, moa = moa)
        for (rep in 2:10) {
            a = cat_outer_fold_pred(res = res, repetition = rep, moa = moa)
            all_rep$data = rbind(all_rep$data, a$data)
        }
        toPlot = generateThreshVsPerfData(all_rep, measures = list( fpr, tpr))
        
        ggplot(toPlot$data, do.call(aes_string,  list(x = "fpr", y = "tpr")) ) +
            geom_path() +
            labs(x = "fpr", y ="tpr") +
            geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5)

    }else{
        dataPlot = c()
        
        for (m in set_MoA){
            all_rep = cat_outer_fold_pred(res = res, repetition = 1, moa = m)
            for (rep in 2:10) {
                a = cat_outer_fold_pred(res = res, repetition = rep, moa = m)
                all_rep$data = rbind(all_rep$data, a$data)
            }
            x = generateThreshVsPerfData(all_rep, measures = list(fpr, tpr))
            x = x$data
            x$MoA = m
            dataPlot = rbind(dataPlot, x)
        }
        toPlot = ggplot(data = dataPlot, do.call(aes_string, list(x = "fpr", y = "tpr"))) +
            geom_path(mapping = aes(color = MoA), size = 2) +
            labs(x = "fpr", y = "tpr") +
            geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5) + 
            scale_color_manual(values = rainbow(length(set_MoA)))
        
        toPlot
    }
}



plot_prec_recall = function(res = result_xgb_5pc_4model_test2, moa = "dna"){
    
    set_MoA = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    
    if(!moa %in% c(set_MoA, "all")){
        print("Invalid Mode of action")
        return(-1)
    }
    if(moa != "all"){
        all_rep = cat_outer_fold_pred(res = res, repetition = 1, moa = moa)
        for (rep in 2:10) {
            a = cat_outer_fold_pred(res = res, repetition = rep, moa = moa)
            all_rep$data = rbind(all_rep$data, a$data)
        }
        toPlot = generateThreshVsPerfData(all_rep, measures = list( ppv, tpr))
        
        ggplot(toPlot$data, do.call(aes_string,  list(x = "tpr", y = "ppv")) ) +
            geom_path() +
            labs(x = "tpr (Recall)", y ="ppv (Precision)")
    }else{
        dataPlot = c()
        
        for (m in set_MoA){
            all_rep = cat_outer_fold_pred(res = res, repetition = 1, moa = m)
            for (rep in 2:10) {
                a = cat_outer_fold_pred(res = res, repetition = rep, moa = m)
                all_rep$data = rbind(all_rep$data, a$data)
            }
            x = generateThreshVsPerfData(all_rep, measures = list(ppv, tpr))
            x = x$data
            x$MoA = m
            dataPlot = rbind(dataPlot, x)
        }
        toPlot = ggplot(data = dataPlot, do.call(aes_string, list(x = "ppv", y = "tpr"))) +
            geom_path(mapping = aes(color = MoA), size = 2) +
            labs(x = "tpr (Recall)", y ="ppv (Precision)") +
            scale_color_manual(values = rainbow(length(set_MoA)))
        
        toPlot
    }
}




#as expected only a few features are really important, most of the time, features convey nearly no information
#Might be interresting to just plot the top X (5,10) features and do a 4 class Venn Diagramm
plot_feat_4model= function(res = result_xgb_5pc_4model_test2, moa = "dna"){
    n_rep = length(res)
    n_folds = length(res[[1]])
   
    all_rep_res = list()
    
    for (i in 1:n_rep){
        rep_res = getFeatureImportance(res[[i]][[1]][[paste0("model_", moa)]])$res
        for (j in 2:n_folds) {
            
            rep_res = rep_res + getFeatureImportance(res[[i]][[j]][[paste0("model_", moa)]])$res
            
        }
        all_rep_res[[i]] = rep_res
    }
    
    aa = do.call(rbind.data.frame, all_rep_res)
    aa = t(aa)
    
    colPal = colorRampPalette(c("white","yellow", "red"), space="rgb") 
    par(oma = c(0,2, 2,1))
    heatmap(aa, scale="column", col = colPal(100) )
    return(sort(apply(aa, 1, sum)))
}

plot_feat_4model()
