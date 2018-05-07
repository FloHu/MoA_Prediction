
#Concatenate prediction of test set from each fold of the outer CV
#Get prediction object for one model of 1vs ALL
#Each individual is predected once
cat_outer_fold_pred = function(res , repetition = 1, moa = "dna"){
    
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
plot_mean_perf = function(res , meas = auc, moa = "dna", ...){
    
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
            pred = cat_outer_fold_pred(res = res,repetition = rep, moa = moa)
            toPlot = c(toPlot, performance(pred, measures = meas))
        } 
        boxplot(toPlot, lwd = 1.5, outline = F, 
                main = paste0("Measure comparison across Repetition :", as.character(meas$id) ) ,
                sub = paste0("MoA : ", moa), ...)
        stripchart(toPlot, method = "jitter", pch = 21, col = "black", bg = rainbow(10), cex = 2, vertical = T, add = T, ...)
    }

}


plot_ROC_allRep = function(res, moa = "dna", plotAllRep = T){
    set_MoA = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    if(!moa %in% c(set_MoA, "all")){
        print("Invalid Mode of action")
        return(-1)
    }
    
    toPlot = ggplot()
    
    if(moa != "all"){
        allData = c()
        
        for (r in 1:10){
            all_test_set = cat_outer_fold_pred(res = res, repetition = r, moa = moa)
            allData = rbind(allData, all_test_set$data)
            if(plotAllRep){
                x = generateThreshVsPerfData(all_test_set, measures = list(fpr, tpr))
                x = x$data
                toPlot = toPlot + geom_path(data = x, mapping = do.call(aes_string, list(x = "fpr", y = "tpr")), size = 0.5, color = "grey", alpha = 0.6)
            }
        }
        all_test_set$data = allData
        plotData = generateThreshVsPerfData(all_test_set, measures = list(fpr, tpr))
        
        toPlot = toPlot + geom_path(data = plotData$data, mapping = do.call(aes_string, list(x = "fpr", y = "tpr")), size = 2) +
                annotate("text", size = 8, x = 0.8, y = 0.1, label = paste0("AUC : ", round(performance(all_test_set, auc), digits = 3)))
    }else{
        colMoa = rainbow(length(set_MoA))
        plotDataMoa = list()
        
        for (m in set_MoA){
            allData = c()
            
            for (r in 1:10){
                all_test_set = cat_outer_fold_pred(res = res, repetition = r, moa = m)
                allData = rbind(allData, all_test_set$data)
            }
            all_test_set$data = allData
            l = paste0(m, " AUC : ", round(performance(all_test_set, auc), digits = 3))
            plotDataMoa[[l]] = all_test_set
           
        }
        plotDataMoa = generateThreshVsPerfData(plotDataMoa, measures = list(fpr, tpr))
        toPlot = plotROCCurves(plotDataMoa) + scale_color_manual(values=colMoa)
        toPlot$layers[[1]]$aes_params$size = 2
        toPlot$labels$colour = "Mode of action"

    }
    toPlot = toPlot +
        geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5) +
        labs(x = "False positive rate", y = "True positive rate", title = paste0(deparse(substitute(res)), " - ", moa )) +
        theme_bw(base_size=15) 
    toPlot
}



#Might be interesting to code a general function taking a list of results 
compare_ROC_2models = function(res1, res2, moa = "dna"){
    set_MoA = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    
    if(!moa %in% c(set_MoA, "all")){
        print("Invalid Mode of action")
        return(-1)
    }
    if(moa != "all"){
        dataPlot = c()
        all_rep_res1 = cat_outer_fold_pred(res = res1, repetition = 1, moa = moa)
        for (rep in 2:10) {
            a = cat_outer_fold_pred(res = res1, repetition = rep, moa = moa)
            all_rep_res1$data = rbind(all_rep_res1$data, a$data)
        }
        all_rep_res2 = cat_outer_fold_pred(res = res2, repetition = 1, moa = moa)
        for (rep in 2:10) {
            a = cat_outer_fold_pred(res = res2, repetition = rep, moa = moa)
            all_rep_res2$data = rbind(all_rep_res2$data, a$data)
        }
        
        
        x1 = generateThreshVsPerfData(all_rep_res1, measures = list( fpr, tpr))
        x1 = x1$data
        x1$Model = deparse(substitute(res1))
        x2 = generateThreshVsPerfData(all_rep_res2, measures = list( fpr, tpr))
        x2 = x2$data
        x2$Model = deparse(substitute(res2))
        dataPlot = rbind(x1, x2)
        
        ggplot(dataPlot, do.call(aes_string,  list(x = "fpr", y = "tpr")) ) +
            geom_path(mapping = aes(color = Model), size = 2) +
            labs(x = "fpr", y ="tpr", title = paste0(deparse(substitute(res)), " - ", moa)) +
            geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5) +
            scale_color_manual(values = rainbow(length(set_MoA))) +
            theme_bw()
        
    }else{
        library(grid)
        library(gridExtra)
        grid.arrange(plot_ROC_allRep(res = res1, moa = "all"), plot_ROC_allRep(res = res2, moa = "all"), nrow = 1)
    }
    
}



plot_prec_recall = function(res , moa = "dna", plotAllRep = T){
    set_MoA = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    if(!moa %in% c(set_MoA, "all")){
        print("Invalid Mode of action")
        return(-1)
    }
    toPlot = ggplot()
    
    if(moa != "all"){
        allData = c()
        for (r in 1:10){
            all_test_set = cat_outer_fold_pred(res = res, repetition = r, moa = moa)
            allData = rbind(allData, all_test_set$data)
            if(plotAllRep){
                x = generateThreshVsPerfData(all_test_set, measures = list(tpr, ppv))
                x = x$data
                toPlot = toPlot + geom_path(data = x, mapping = do.call(aes_string, list(x = "tpr", y = "ppv")), size = 0.5, color = "grey", alpha = 0.6)
            }
        }
        all_test_set$data = allData
        plotData = generateThreshVsPerfData(all_test_set, measures = list(tpr, ppv))
        
        toPlot = toPlot +
            geom_path(data = plotData$data, mapping = do.call(aes_string, list(x = "tpr", y = "ppv")), size = 2)
        
    }else{
        colMoa = rainbow(length(set_MoA))
        plotDataMoa = list()
        
        for (m in set_MoA){
            allData = c()
            
            for (r in 1:10){
                all_test_set = cat_outer_fold_pred(res = res, repetition = r, moa = m)
                allData = rbind(allData, all_test_set$data)
            }
            all_test_set$data = allData
            plotDataMoa[[m]] = all_test_set
            
        }
        plotDataMoa = generateThreshVsPerfData(plotDataMoa, measures = list(tpr, ppv))
        toPlot = plotROCCurves(plotDataMoa, diagonal = F ) + scale_color_manual(values=colMoa)
        toPlot$layers[[1]]$aes_params$size = 2
        toPlot$labels$colour = "Mode of action"
        
    }
    toPlot = toPlot +  labs(x = "tpr (Recall)", y ="ppv (Precision)", title = paste0(deparse(substitute(res)), " - ", moa )) + theme_bw()
    toPlot
}




#as expected only a few features are really important, most of the time, features convey nearly no information
#Might be interresting to just plot the top X (5,10) features and do a 4 class Venn Diagramm
plot_feat_4model= function(res , moa = "dna"){
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

    


plot_predProb_moa = function(res, rep, moa = "dna"){
    set_MoA = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    if(!moa %in% c(set_MoA, "all")){
        print("Invalid Mode of action")
        return(-1)
    }
    
    dt = cat_outer_fold_pred(res =res, repetition = rep, moa = moa)
    dt = dt$data
    
    dt_moa = filter(dt, truth == moa)
    dt_not_moa = filter(dt, truth == paste0("not_", moa))
    
    toPlot = list(dt_moa[, paste0("prob.", moa)], dt_moa[, paste0("prob.not_", moa)], dt_not_moa[, paste0("prob.", moa)], dt_not_moa[, paste0("prob.not_", moa)])
    boxplot(toPlot, at = c(1,2,4,5), main = "Truth", col = rainbow(2), lwd = 2, xaxt = "n")
    axis(side = 1, at = c(1,2,4,5), labels = c(paste0("prob ", moa), paste0("not_", moa), paste0("prob ", moa), paste0("not_", moa)))
    axis(side = 3, at = c(1.5,4.5), labels = c(moa, paste0("not_", moa)))
    
}





# ==============================================================================



plot_ROC_optThres = function(res, moa = "dna"){
    
    if(!require(ggrepel)){
        install.packages("ggrepel")
    }
    library(ggrepel)
    library(grid)
    library(gridExtra)
    
    #Concatenate all Data from all repetition/outerFolds
    allData = c()
    for (r in 1:10){
        all_test_set = cat_outer_fold_pred(res = res, repetition = r, moa = moa)
        allData = rbind(allData, all_test_set$data)
    }
    all_test_set$data = allData
    
    
    mccCurve = generateThreshVsPerfData(all_test_set, measures = mcc)
    #Best Threshold
    bestThres = mccCurve$data[which.max(mccCurve$data$mcc),"threshold"]
    
    #Custom Threshold
    all_test_set = setThreshold(all_test_set, threshold = bestThres)
    confMat = calculateConfusionMatrix(all_test_set)
    TP = confMat$result[1,1]
    FP = confMat$result[2,1]
    FN = confMat$result[1,2]
    TN = confMat$result[2,2]
    
    
    plotData = generateThreshVsPerfData(all_test_set, measures = list(fpr, tpr))
    toPlot = ggplot() + geom_path(data = plotData$data, mapping = do.call(aes_string, list(x = "fpr", y = "tpr")), size = 2) +
        annotate("text", size = 8, x = 0.8, y = 0.1, label = paste0("AUC : ", round(performance(all_test_set, auc), digits = 3)))
    
    
    
    bestPoint = plotData$data[which(round(plotData$data$threshold, digits = 2) == round(bestThres, digits = 2)) , ]
    
    toPlot = toPlot + annotate("point", x = bestPoint$fpr, y = bestPoint$tpr, colour = "red", cex = 5) + 
        geom_label(aes(label = paste0("Threshold = ", round(bestThres, digits = 2) ,"\nTP : ", TP, " FP : ", FP, "\nFN : ", FN, " TN : ", TN),
                            size = 3.5, x = bestPoint$fpr+0.15, y = bestPoint$tpr - 0.05 ), show.legend = FALSE) +
        geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5) +
        labs(x = "False positive rate", y = "True positive rate", title = paste0(deparse(substitute(res)), " - ", moa )) 
        
    

    grid.arrange(plotThreshVsPerf(mccCurve), toPlot , nrow = 1)
}
