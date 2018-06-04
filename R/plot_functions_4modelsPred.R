# ==============================================================================
#
#              FUNCTIONS GENERATING PLOTS FROM RESULTS OBJECTS
#
# ==============================================================================

#AIM : Concatenate prediction of test set from each fold of the outer CV,
#       so each individual is predected once
#INPUT : Result object, repetition number and moa name
#OUTPUT : Prediction object countaining all individuals predicted once
#           with a model that has not been train on them
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

# ==============================================================================

#AIM : Generate a boxplot of a performance measure across all repetition
#       Can be used for 1 MoA or for comparing all of them
#INPUT : Result object, any measure from the MLR package and a MoA name or "all"
#OUTPUT : Some classic boxplot
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
        boxplot(t(toPlot), lwd = 1.5, outline = F, ylab = meas$name,
                main = paste0("Measure comparison across MoA :", as.character(meas$id) ),
                sub = paste0("result object : ", deparse(substitute(res))),
                col = rainbow(length(set_MoA)),...)
    }else{
        toPlot = c()
        for(rep in 1:10){
            pred = cat_outer_fold_pred(res = res,repetition = rep, moa = moa)
            toPlot = c(toPlot, performance(pred, measures = meas))
        }
        boxplot(toPlot, lwd = 1.5, outline = F, ylab = meas$name,
                main = paste0("Measure comparison across Repetition : ", as.character(meas$id), " for ", moa, " MoA" ) ,
                sub = paste0("result object : ", deparse(substitute(res))), ...)
        stripchart(toPlot, method = "jitter", pch = 21, col = "black", bg = "orange", cex = 2, vertical = T, add = T,lwd = 1.5, ...)
    }
}

# ==============================================================================

#AIM : Overall results of a run, plot a ROC curve based on the concatenation of
#       all outer folds of all repetitions
#INPUT : Results object, MoA name or "all"
#OUTPUT : ROC curve(s) with AUC as a legend
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
            l = paste0(m, " AUC : ", round(performance(pred = all_test_set, measures = auc), digits = 3))
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

# ==============================================================================

#AIM : Compare thr ROC curves beetween two results object
#INPUT : 2 result objects res1 and res2, 1 MoA name or "all"
#OUTPUT : One plot with one ROC curve per Model if 1 MoA
#        A splited plot with the ouput of plot_ROC_allRep() for each model
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
            labs(x = "False positive rate", y = "True positive rate", title = paste0(moa, " - Comparison")) +
            geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5) +
            scale_color_manual(values = rainbow(length(set_MoA))) +
            theme_bw()

    }else{
        library(grid)
        library(gridExtra)
        plot1 = plot_ROC_allRep(res = res1, moa = "all") +
            theme(legend.position = "bottom", legend.box = "vertical") +
            guides(col = guide_legend(title.position = "top", nrow = 2)) +
            ggtitle(deparse(substitute(res1)))
        plot2 = plot_ROC_allRep(res = res2, moa = "all") +
            theme(legend.position = "bottom", legend.box = "vertical") +
            guides(col = guide_legend(title.position = "top", nrow = 2)) +
            ggtitle(deparse(substitute(res2)))

        grid.arrange(plot1, plot2, nrow = 1)
    }
}


#Same thing as above with any amount of models
#Work only for one MoA
compare_ROC_models = function(moa = "dna", ...){

    if(!moa %in% c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
        print("Invalid Mode of action")
        return(-1)
    }
    
    res_objs = list(...)
    names(res_objs) = sapply(substitute(list(...))[-1], deparse)

    if(length(res_objs) == 0){
        print("No models given in arguments")
        return(-1)
    }
    
    dataPlot = c()
    for (i in 1:length(res_objs)) {
        all_rep_res = cat_outer_fold_pred(res = res_objs[[i]], repetition = 1, moa = moa)
        for (rep in 2:10) {
            a = cat_outer_fold_pred(res = res_objs[[i]], repetition = rep, moa = moa)
            all_rep_res$data = rbind(all_rep_res$data, a$data)
        }
        x1 = generateThreshVsPerfData(all_rep_res, measures = list( fpr, tpr))
        x1 = x1$data
        x1$Model = paste0(names(res_objs)[i], " AUC : ", round(performance(all_rep_res, auc), digits = 3))
        dataPlot = rbind(dataPlot, x1)
    }
    dataPlot$Model = factor(x = dataPlot$Model, levels = unique(x = dataPlot$Model))

    ggplot(dataPlot, do.call(aes_string,  list(x = "fpr", y = "tpr")) ) +
        geom_path(mapping = aes(color = Model), size = 2) +
        labs(x = "False positive rate", y = "True positive rate", title = paste0(moa, " - Comparison across models")) +
        geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5) +
        scale_color_manual(values = rainbow(length(res_objs))) +
        theme_bw()
}

# ==============================================================================

#AIM : Results of a run, plot a Precision-Recall curve based on the concatenation of
#       all outer folds of all repetitions
#INPUT : Results object, MoA name or "all"
#OUTPUT : Precision-Recall curve
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

# ==============================================================================

#AIM : Return Features Importance across all outer fold of all repetition
#INPUT : Result object and 1 MoA name
#OUTPUT : Heatmap plot of feature importance + return names vector of the sumed
#           features importances
#NB : Does not work with wilcoxon based feature selection because all folds don't
#       have the same number of features
plot_feat_4model= function(res , moa = "dna", noPlot = F){
    if(!moa %in% c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
        print("Invalid Mode of action")
        return(-1)
    }
    n_rep = length(res)
    n_folds = length(res[[1]])

    all_rep_res = list()
    for (i in 1:n_rep){
        rep_res = getFeatureImportance(res[[i]][[1]][[paste0("model_", moa)]])$res
        for (j in 2:n_folds){
            rep_res = rep_res + getFeatureImportance(res[[i]][[j]][[paste0("model_", moa)]])$res
        }
        all_rep_res[[i]] = rep_res
    }

    aa = do.call(rbind.data.frame, all_rep_res)
    aa = t(aa)

    if(!noPlot){
        colPal = colorRampPalette(c("white","yellow", "red"), space="rgb")
        par(oma = c(0,2, 2,1))
        heatmap(aa, scale="column", col = colPal(100), main = paste0("Features importance : ", moa), sub = paste0("dataset", deparse(substitute(res))),
      cexRow = 0.3)
      plot <- recordPlot()
    } else {
      plot <- NA
    }
    return(list(features = sort(apply(aa, 1, sum)),
                plot = plot))
}

# ==============================================================================

#AIM : Plot the distribution of prediction probabilities for drugs labelled MoA and Not_MoA
#INPUT : Result object, 1 MoA name, one repetition number
#OUTPUT : Boxplot of probabilities distribution by class
plot_predProb_moa = function(res, rep = 1, moa = "dna", dt_matrix = the_matrix_allDrugs){
    if(!moa %in% c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
        print("Invalid Mode of action")
        return(-1)
    }
    color_moa = c(rainbow(4), rep("black", 3))
    names(color_moa) = c("dna", "cell_wall", "membrane_stress", "protein_synthesis", "protein_qc", "oxidative_stress", "pmf")

    dt = cat_outer_fold_pred(res =res, repetition = rep, moa = moa)
    dt = dt$data

    dt_moa = filter(dt, truth == moa)
    dt_not_moa = filter(dt, truth == paste0("not_", moa))

    toPlot = list(dt_moa[, paste0("prob.", moa)], dt_moa[, paste0("prob.not_", moa)], dt_not_moa[, paste0("prob.", moa)], dt_not_moa[, paste0("prob.not_", moa)])
    boxplot(toPlot, at = c(1,2,4,5), main = "Truth", col = rainbow(2), lwd = 2, xaxt = "n", outline = F,
            sub = paste0("Dataset : ", deparse(substitute(res)), " - Rep ", rep))

    stripchart(at = 4, vertical = T, dt_not_moa[, paste0("prob.", moa)], pch = 21, bg = color_moa[the_matrix_allDrugs[dt_not_moa$id, "process_broad"]] , cex = 1.2, method = "jitter", add =T)
    stripchart(at = 5, vertical = T, dt_not_moa[, paste0("prob.not_", moa)], pch = 21, bg = color_moa[the_matrix_allDrugs[dt_not_moa$id, "process_broad"]] , cex = 1.2, method = "jitter", add =T)

    axis(side = 1, at = c(1,2,4,5), labels = c(paste0("prob ", moa), paste0("not_", moa), paste0("prob ", moa), paste0("not_", moa)))
    axis(side = 3, at = c(1.5,4.5), labels = c(moa, paste0("not_", moa)))
}

# ==============================================================================

#AIM : Plot a ROC curve for one MoA as well as the MCC depending on the threshold
#       For the optimal threshold (highest MCC), display a confusion matrix on the ROC curve
#INPUT :  A result object and a MoA
#OUTPUT : Splitted double plot of the 2 curves and return a confusionMatrix object
plot_ROC_optThres = function(res, moa = "dna", customTitle = paste0(deparse(substitute(res)), " - ", moa )){
    
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
    
    defaultTPR = performance(all_test_set, measures = tpr)
    defaultFPR = performance(all_test_set, measures = fpr)
    defaultConfMat = calculateConfusionMatrix(all_test_set)
    TP2 = defaultConfMat$result[1,1]
    FP2 = defaultConfMat$result[2,1]
    FN2 = defaultConfMat$result[1,2]
    TN2 = defaultConfMat$result[2,2]
    
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
        annotate("text", size = 10, x = 0.8, y = 0.1, label = paste0("AUC : ", round(performance(all_test_set, auc), digits = 3))) +
        theme(text = element_text(size = 16))
    
    bestPoint = plotData$data[which(round(plotData$data$threshold, digits = 2) == round(bestThres, digits = 2)) , ]
    
    toPlot = toPlot + annotate("point", x = bestPoint$fpr, y = bestPoint$tpr, colour = "red", cex = 5) +
        annotate("point", x = defaultFPR, y = defaultTPR, colour = "blue", cex = 5) +
        geom_label(aes(label = paste0("Threshold = 0.5 (blue)\nTP : ", TP2, " FP : ", FP2, "\nFN : ", FN2, " TN : ", TN2),
                       size = 3.5, x = 0.85, y = 0.4 ), show.legend = FALSE) +
        geom_label(aes(label = paste0("Threshold = ", round(bestThres, digits = 2) ," (red)\nTP : ", TP, " FP : ", FP, "\nFN : ", FN, " TN : ", TN),
                       size = 3.5, x = 0.85, y = 0.2 ), show.legend = FALSE) +
        geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5) +
        labs(x = "False positive rate", y = "True positive rate", title = customTitle) +
        theme(text = element_text(size = 16))
    
    toPlot <- grid.arrange(plotThreshVsPerf(mccCurve) + geom_vline(xintercept=bestThres, linetype = "dotted"), toPlot , nrow = 1)
    return(list(confMat = confMat, plot = toPlot))
}


# ==============================================================================

#AIM : Display distribution of prediction probabilities by MoA for one drug in a result object
#INPUT : Result object, drugname (typaslab)
#OUTPUT : Boxplot of probabilities across all repetitions
distrib_drug_prob = function(res, drug = "A22", dt_matrix = the_matrix_allDrugs){

    idDrug =  which(dt_matrix$drugname_typaslab == drug)
    drugMoa = filter(the_matrix_allDrugs, drugname_typaslab == drug) %>% select(process_broad)

    prob_allMoa = list()

    for (m in c("cell_wall", "dna", "membrane_stress", "protein_synthesis")) {
        allData = c()
        for (r in 1:10){
            all_test_set = cat_outer_fold_pred(res = res, repetition = r, moa = m)
            allData = rbind(allData, all_test_set$data)
        }
        all_test_set$data = allData
        prob_allMoa[[paste0("prob_", m)]] = filter(all_test_set$data, id == idDrug ) %>% select(paste0("prob.", m)) %>% t(.)
    }
    boxplot(prob_allMoa, lwd = 2, col = rainbow(length(prob_allMoa)),
            main = paste0("Prediction Probabilities - ", drug, " - Truth = ", drugMoa),
            sub = paste0("Dataset : ", deparse(substitute(res))))
    stripchart(prob_allMoa, vertical = T, add = T, method = "jitter", pch = 21, bg = "grey", cex = 1.5)
}
