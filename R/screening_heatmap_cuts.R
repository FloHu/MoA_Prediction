



recall_from_table = function(tab, moa, mat){
    return(tab[moa] / sum(mat$process_broad == moa))
}
# ==============================================================================


ppv_from_table = function(tab, moa){
    return(tab[moa] / sum(tab))
}
# ==============================================================================


screening_heatmap_cuts = function(res_obj_file, model_type, moa, mat_from_container){
    
    library(ComplexHeatmap)
    library("circlize")
    
    colMap = rainbow(4)
    names(colMap) = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    res_obj = readRDS(file = res_obj_file)
    res_cuts_perf = list()
    
    if(model_type == "lasso"){
        seq_thres = seq(from = 1, to = 0.5, by = -0.1)
    }else{
        seq_thres = seq(from = 1, to = 0.8, by = -0.0125)
    }
    
    for (thres in seq_thres) {
        cat("Threshold =", thres, "\n")
        featImp = plot_top_feat_importance(resObj = res_obj, moa = moa, main = deparse(substitute(resObj)) , thres = thres, return_obj = T, model_type = model_type)  
        
        cluster_matrix = mat_from_container %>% select(drugname_typaslab, process_broad, names(featImp))
        cluster_matrix = as.data.frame(cluster_matrix)
        rownames(cluster_matrix) = paste0(cluster_matrix$process_broad, seq(1, nrow(cluster_matrix)))
        
        annot = data.frame(MoA = cluster_matrix[ ,"process_broad"])
        rownames(annot) = rownames(cluster_matrix)
        
        cluster_matrix = cluster_matrix %>% select(-drugname_typaslab)
        
        if(ncol(cluster_matrix) <= 7){
            seq_geneCut = 1:(ncol(cluster_matrix)-1)
        }else{
            seq_geneCut = 1:6
        }
        
        for(geneCut in seq_geneCut){
            cat("Gene Cut =", geneCut, "\n")
            # First cut : gene groups (columns)
            genes_grp = cutree(hclust(dist(t(cluster_matrix %>% select(-"process_broad")))), k = geneCut)
            mats_subGene = list()
            # Defining sub matrices based on gene groups
            for (i in 1:geneCut){
                mats_subGene[[i]] = cluster_matrix[ ,c(names(genes_grp[genes_grp == i]), "process_broad")]
            }
            for (i in 1:length(mats_subGene)){
                # For every submatrix of gene, we cut drugs groups (rows hclust)
                for (drugsCut in 2:6){
                    cat("Drugs Cut =", drugsCut, "\n")
                    grp = cutree(hclust(dist(mats_subGene[[i]] %>% select(-"process_broad"))), k = drugsCut)
                    for (g in 1:drugsCut){
                        # For every group defined, assess the quality of it in term of recall and ppv
                        t = table(mats_subGene[[i]][names(grp[grp == g]), "process_broad"])
                        res_cuts_perf = c(res_cuts_perf, list(threshold = thres, nGeneCut = geneCut, nDrugsCut = drugsCut, 
                                          recall = recall_from_table(tab = t, moa = moa, mat = cluster_matrix),
                                          ppv = ppv_from_table(tab = t, moa = moa), tab_value = t[moa], nbGene = ncol(mats_subGene[[i]]))
                                        )
                    }
                }
            }
        }
   
    }
    
    res_cuts_perf = data.frame(matrix(unlist(res_cuts_perf), ncol = 7, byrow = T))
    colnames(res_cuts_perf) = c("threshold", "nGeneCut", "nDrugsCut", "recall", "ppv", "tab", "nbGene")
    
    p1 = ggplot(data = res_cuts_perf, aes(ppv, recall)) + geom_point(size = 2, alpha = 0.5) + theme_bw()
    p2 = ggplot(data = res_cuts_perf, aes(ppv, recall)) + geom_jitter(size = 2, alpha = 0.5, width = 0.01, height = 0.01) + theme_bw()
    grid.arrange(p1, p2)
   
    #Compute distance to optimal point : 1 recall and 1 PPV = full and pure cluster of targeted MoA
    res_cuts_perf$distToOptPoint = sqrt((1-res_cuts_perf$recall)**2 + (1-res_cuts_perf$ppv)**2)
    #res_cuts_perf = na.omit(res_cuts_perf)
    

    return(res_cuts_perf)
}
# ==============================================================================



