
plot_heatmap <- function(m  , plotFile = F, name_size = 0.3 , names_order = NULL, scaleMethod = "column", clust = "row") {
    # plots a heatmap for all drugs with only the features provided in the feats_to_keep argument
    cluster_matrix <- m
    moa_to_colour <- c(cell_wall = "#FF0000FF",
                       dna = "#80FF00FF",
                       protein_synthesis = "#8000FFFF",
                       membrane_stress = "#00FFFFFF",
                       pmf = "#bababa", protein_qc = "#bababa", oxidative_stress = "#bababa"
    )
    names_order = rev(names_order)

    if(!is.null(names_order)){
        row_ord = unlist(lapply(names_order, function(x){which(cluster_matrix$drugname_typaslab== x)}))
    }else{
        row_ord = TRUE
    }
    drugs_names = cluster_matrix$drugname_typaslab[row_ord]
    cluster_matrix = cluster_matrix[ row_ord, ]
    
    process_lut <- cluster_matrix$process_broad
    
    rowcols <- moa_to_colour[process_lut]
   
    
    cluster_matrix$drugname_typaslab <- NULL
    cluster_matrix$process_broad <- NULL
    cluster_matrix <- as.matrix(cluster_matrix)
    
    if(plotFile){ pdf("./plots/heatmap.pdf", width = 20, height = 20) }
    
    if(is.null(clust)){
        row_ord = FALSE
        dend = "none"
    }else{
        dend = clust
    }
    
    heatmap.2(cluster_matrix,
              scale = scaleMethod,
              trace = "none",
              Rowv = row_ord, 
              Colv = TRUE,
              breaks = c(seq(-5,5,length=20)), # use this i/o zlim
              RowSideColors = rowcols,
              cexRow = name_size,
              labRow = drugs_names,
              col = colorRampPalette(c("red", "white", "blue")),
              margins = c(12, 9),
              dendrogram = dend) # only draw row dendrogram
    legend("top",
           legend = names(moa_to_colour),
           col = moa_to_colour,
           lty = 1,
           lwd = 10)
    
    if(plotFile){ dev.off() }
}

# =========================================================================================================================================

plot_top_feat_importance = function(resObj, moa = "dna", thres = 0, return_obj = FALSE, model_type = "lasso", ...){
    
    match.arg(arg = moa, choices = c("cell_wall", "dna", "membrane_stress", "protein_synthesis"))
    match.arg(arg = model_type, choices = c("lasso", "tree"))
    
    nrep = length(resObj)
    nfold = length(resObj[[1]])
    #List of coefficient values across all models
    feat_cat = list()
    
    for(rep in 1:nrep){
        for(fold in 1:nfold){
            
            if(model_type == "lasso"){
                model = resObj[[rep]][[nfold]][[paste0("model_", moa)]]
                lambda_for_pred = model$learner$par.vals$s
                # Get the index of the lambda used in model building that is the closest to s, the lambda used for testing
                closest_lambda_index = which.min(abs(lambda_for_pred - model$learner.model$lambda))
                coeffs = as.matrix(model$learner.model$beta[, closest_lambda_index])
                coeffs = coeffs[coeffs !=0,  ]
            }else if(model_type == "tree"){
                coeffs = getFeatureImportance(resObj[[rep]][[fold]][[paste0("model_", moa)]])$res
                coeffs = coeffs[which(coeffs !=0)]
            }
            
            for(n in names(coeffs)){
                feat_cat[[n]] = unlist(c(feat_cat[[n]], coeffs[n]))
            }
        }
    }
    
    #eff contains how many times features appear in a model
    eff = sort(unlist(lapply(feat_cat, length)))
    feat_cat = feat_cat[names(eff)]
    
    # filter
    feat_cat = feat_cat[(eff/80) >= thres]
    eff = eff[(eff/80) >= thres]
    
    colMap = rainbow(length(unique(eff)))
    names(colMap) = unique(eff)
    par(mar = c(5,10,4,2))
    boxplot(feat_cat, horizontal = T, las = 2, lwd = 1.5, col = colMap[as.character(eff)], ...)
    
    if(model_type == "lasso"){
       location = "bottomleft"
    }else if(model_type == "tree"){
        location = "bottomright"
    }
    
    legend(x = location,  legend = paste0(as.numeric(names(colMap))/80 * 100, " %"), 
           fill = colMap, title = "Present in % models", cex = 0.8)
    abline(v = 0)
    if(return_obj){
        return(feat_cat)
    }
}

# =========================================================================================================================================

model_analysis = function(res_obj , matrix_container_line, matrix_container_line_noChemFeat, moa,
                          pdf_filename = "plot.pdf", feat_imp_thres = 0.6, model_type = "lasso", topFeatThres = 0.9){
    
    nb_models = length(res_obj) * length(res_obj[[1]])
    mat = matrix_container_line$drug_feature_matrices[[1]]
    colMap = rainbow(4)
    names(colMap) = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    
    # PDF file with all plot from the analysis
    pdf(file = pdf_filename, width = 14, height = 14, onefile = T)
    
    #First page : ROC curve of the model for the targeted MoA
    par(mfrow = c(1,1))
    grid.arrange(plot_ROC_allRep(res_obj, moa = moa, plotAllRep = T))
    grid.arrange(plot_prec_recall(res_obj, moa = moa, plotAllRep = T))
    
    featImp = plot_top_feat_importance(resObj = res_obj, moa = moa, main= deparse(substitute(resObj)) , thres = feat_imp_thres, return_obj = T, model_type = model_type)    
    a = unlist(lapply(featImp, length))
    topFeat = names(a[a >= nb_models*topFeatThres])
   
    # Stripcharts of top feat repartition
    par(mfrow = c(2,2))
    for(f in topFeat){
        stripchart(mat[[f]] ~ mat$process_broad, pch = 21, bg = "orange", cex = 1.2, las = 2, xlab = f)
    }
    
    # Plot of top features 2 by 2
    for(f1 in 1:(length(topFeat)-1)){
        for(f2 in 2:length(topFeat)){
            if(f1 != f2){
                plot(x = mat[[topFeat[f1]]], y = mat[[topFeat[f2]]], pch = 21, bg = colMap[mat$process_broad], cex = 1.5, las = 2, xlab = topFeat[f1], ylab = topFeat[f2])
            }
        }
    }
    
    par(mfrow = c(1,2))
    mat = mat %>% select(drugname_typaslab, process_broad, names(featImp))
    names_order = drugs_pred_prob_from_container(pred_data = matrix_container_line$PredData[[1]], moa = moa,
                                                main = paste0(matrix_container_line$hyperparam_grid_name,
                                                              matrix_container_line$feat_preselect,
                                                              as.character(matrix_container_line$chemical_feats)), sep = "_" )
    prob_shift_plot(pred_data = matrix_container_line$PredData[[1]], pred_data_noChem = matrix_container_line_noChemFeat$PredData[[1]], moa = moa, main = "Comparison with or without ChemFeat" )
    
    
    if(nrow(mat) > 100){
        plot_heatmap(m = mat, names_order = names_order)
        plot_heatmap(m = mat, names_order = names_order, scaleMethod = "none")
        plot_heatmap(m = mat, names_order = names_order, clust = "column")
    }else{
        plot_heatmap(m = mat, names_order = names_order, name_size = 0.9)
        plot_heatmap(m = mat, names_order = names_order, scaleMethod = "none", name_size = 0.9)
        plot_heatmap(m = mat, names_order = names_order, clust = "column", name_size = 0.9)
    }
   
    
    dev.off()
}

