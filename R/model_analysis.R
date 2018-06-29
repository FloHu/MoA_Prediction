
plot_heatmap <- function(m  , plotFile = F, name_size = 0.3 , names_order = NULL) {
    # plots a heatmap for all drugs with only the features provided in the feats_to_keep argument
    cluster_matrix <- m
    moa_to_colour <- c(cell_wall = "#FF0000FF",
                       dna = "#80FF00FF",
                       protein_synthesis = "#8000FFFF",
                       membrane_stress = "#00FFFFFF",
                       pmf = "#bababa", protein_qc = "#bababa", oxidative_stress = "#bababa"
    )
    
    process_lut <- cluster_matrix$process_broad
    
    rowcols <- moa_to_colour[process_lut]
    drugs_names = cluster_matrix$drugname_typaslab
    if(!is.null(names_order)){
        row_ord = unlist(lapply(names_order, function(x){which(cluster_matrix$drugname_typaslab== x)}))
    }else{
        row_ord = TRUE
    }
    
    cluster_matrix$drugname_typaslab <- NULL
    
    cluster_matrix$process_broad <- NULL
    cluster_matrix <- as.matrix(cluster_matrix)
    
    if(plotFile){ pdf("./plots/heatmap.pdf", width = 20, height = 20) }
    
    heatmap.2(cluster_matrix,
              scale = "column",
              trace = "none",
              Rowv = row_ord, 
              breaks = c(seq(-5,5,length=20)), # use this i/o zlim
              RowSideColors = rowcols,
              cexRow = name_size,
              labRow = drugs_names,
              col = colorRampPalette(c("red", "white", "blue")),
              margins = c(12, 9),
              dendrogram = "row") # only draw row dendrogram
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
            }else if(model == "tree"){
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
    }else if(model == "tree"){
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

model_analysis = function(res_obj , matrix_container_line , moa, pdf_filename = "plot.pdf", feat_imp_thres = 0.6, model_type = "lasso"){
    
    nb_models = length(res_obj) * length(res_obj[[1]])
    mat = matrix_container_line$drug_feature_matrices[[1]]
    colMap = rainbow(4)
    names(colMap) = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    
    # PDF file with all plot from the analysis
    pdf(file = pdf_filename, width = 14, height = 14 )
    
    #First page : ROC curve of the model for the targeted MoA
    par(mfrow = c(1,1))
    plot_ROC_allRep(res_obj, moa = moa, plotAllRep = T)
    
    featImp = plot_top_feat_importance(resObj = res_obj, moa = moa, main= deparse(substitute(resObj)) , thres = feat_imp_thres, return_obj = T, model_type = model_type)    

    
    a = unlist(lapply(featImp, length))
    topFeat = names(a[a >= nb_models*0.9])
    
    # Stripcharts of top feat repartition
    par(mfrow = c(2,2))
    for(f in topFeat){
        stripchart(mat[[f]] ~ mat$process_broad, pch = 21, bg = "orange", cex = 1.2, las = 2, xlab = f)
    }
    
    # Plot of top features 2 by 2
    for(f1 in topFeat){
        for(f2 in topFeat){
            if(f1 != f2){
                plot(mat[[f1]] ~ mat[[f2]], data = mat, pch = 21, bg = colMap[process_broad], cex = 1.5, las = 2, xlab = f1, ylab = f2)
            }
        }
    }
    
    par(mfrow = c(1,1))
    mat = mat %>% select(drugname_typaslab, process_broad, names(featImp))
    names_order = drugs_pred_prob_from_container(pred_data = matrix_container_line$PredData[[1]], moa = moa )
    plot_heatmap(m = mat, names_order = names_order)
    
    dev.off()
}

