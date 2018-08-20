



recall_from_table = function(tab, moa, mat){
    if(is.na(tab[moa])){
        return(0)
    }else{
        return(tab[moa] / sum(mat$process_broad == moa))
    }
}
# ==============================================================================


ppv_from_table = function(tab, moa){
    if(is.na(tab[moa])){
        return(0)
    }else{
        return(tab[moa] / sum(tab))
    }
}
# ==============================================================================


screening_heatmap_cuts = function(res_obj_file, model_type, moa, mat_from_container, dist_method = "euclidian"){
    
    library(ComplexHeatmap)
    library("circlize")
    
    colMap = rainbow(4)
    names(colMap) = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    res_obj = readRDS(file = res_obj_file)
    res_cuts_perf = list()
    

    seq_thres = seq(from = 1, to = 0.85, by = -0.0125)

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
            if(dist_method == "correlation_based"){
                dist_corr = sqrt(2*(1 - cor(cluster_matrix %>% select(-"process_broad"))))
                genes_grp = cutree(hclust(as.dist(dist_corr)), k = geneCut)
            }else{
                genes_grp = cutree(hclust(dist(t(cluster_matrix %>% select(-"process_broad")))), k = geneCut)
            }
            mats_subGene = list()
            # Defining sub matrices based on gene groups
            for (i in 1:geneCut){
                mats_subGene[[i]] = cluster_matrix[ ,c(names(genes_grp[genes_grp == i]), "process_broad")]
            }
            for (i in 1:length(mats_subGene)){
                if(ncol(mats_subGene[[i]]) <= 3){
                    next
                }
                # For every submatrix of gene, we cut drugs groups (rows hclust)
                for (drugsCut in 2:6){
                    cat("Drugs Cut =", drugsCut, "\n")
                    
                    if(dist_method == "correlation_based"){
                        dist_corr = sqrt(2*(1 - cor(t(mats_subGene[[i]] %>% select(-"process_broad")))))
                        grp = cutree(hclust(as.dist(dist_corr)), k = drugsCut)
                    }else{
                        grp = cutree(hclust(dist(mats_subGene[[i]] %>% select(-"process_broad"))), k = drugsCut)
                    }
                    
                    for (g in 1:drugsCut){
                        # For every group defined, assess the quality of it in term of recall and ppv
                        t = table(mats_subGene[[i]][names(grp[grp == g]), "process_broad"])
                        res_cuts_perf = c(res_cuts_perf, list(threshold = thres, nGeneCut = geneCut, nDrugsCut = drugsCut, 
                                          recall = recall_from_table(tab = t, moa = moa, mat = cluster_matrix),
                                          ppv = ppv_from_table(tab = t, moa = moa), 
                                          tab_value = ifelse(is.na(t[moa]), yes = 0, no = t[moa]),
                                          nbGene = ncol(mats_subGene[[i]]))
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
    
    rm(res_obj)
    gc()
    return(res_cuts_perf)
}


# ==============================================================================
# Lots of loops, testing different way of splitting every result object to find the best feat set
# Might take a while
# Might also crash

# Actually take quite a while (like several hours), don't use it

screening_heatmap_cuts_HARDCORE_VERSION = function( res_file_pattern = "^most2", moa, full_matrix){
    
    library(ComplexHeatmap)
    library("circlize")
    
    colMap = rainbow(4)
    names(colMap) = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    res_cuts_perf = list()
    
    
    res_files = list.files(path = "run_results_from_server/matrix_container_result/", pattern = res_file_pattern)
    
    # FOR EACH RES FILE
    for(file in res_files){
        res_obj = readRDS(file = paste0("run_results_from_server/matrix_container_result/", file) )
        
        # ASSESS MODEL TYPE
        if(grepl(x = file, pattern = "classif.glmnet")){
            model_type = "lasso"
            seq_thres = seq(from = 1, to = 0.6, by = -0.1)
        }else{
            model_type = "tree"
            seq_thres = seq(from = 1, to = 0.8, by = -0.0125)
        }
        
        # SCREEN THRESHOLD VALUES
        for (thres in seq_thres) {
            cat("Threshold =", thres, "\n")
            featImp = plot_top_feat_importance(resObj = res_obj, moa = moa, main = deparse(substitute(resObj)) , thres = thres, return_obj = T, model_type = model_type)  
            
            cluster_matrix = full_matrix %>% select(drugname_typaslab, process_broad, names(featImp))
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
                for (distMethod in c("dist_euclidian", "pearson_corr")) {
                    cat("Gene Cut =", geneCut, "\n")
                    
                    # First cut : gene groups (columns)
                    if(distMethod == "dist_euclidian"){
                        genes_grp = cutree(hclust(dist(t(cluster_matrix %>% select(-"process_broad")))), k = geneCut)
                    }else{
                        corTree = hclust(as.dist(cor(cluster_matrix %>% select(-process_broad))))
                        genes_grp = cutree(corTree, k = geneCut)
                    }
                    
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
                                res_cuts_perf = c(res_cuts_perf, list(file = file,
                                                                    distMethod = distMethod,
                                                                    threshold = thres, nGeneCut = geneCut, nDrugsCut = drugsCut, 
                                                                    recall = recall_from_table(tab = t, moa = moa, mat = cluster_matrix),
                                                                    ppv = ppv_from_table(tab = t, moa = moa), 
                                                                    tab_value = ifelse(is.na(t[moa]), yes = 0, no = t[moa]),
                                                                    nbGene = ncol(mats_subGene[[i]]) -1)
                                )
                            }
                        }
                    } 
                }
  
            }
        }
        #Cleaning for memory optimization
        rm(res_obj)
        gc()
    }
    
    res_cuts_perf = data.frame(matrix(unlist(res_cuts_perf), ncol = 9, byrow = T))
    colnames(res_cuts_perf) = c("file", "distMethod", "threshold", "nGeneCut", "nDrugsCut", "recall", "ppv", "tab", "nbGene")
    
    p1 = ggplot(data = res_cuts_perf, aes(ppv, recall)) + geom_point(size = 2, alpha = 0.3) + theme_bw()
    p2 = ggplot(data = res_cuts_perf, aes(ppv, recall)) + geom_jitter(size = 2, alpha = 0.5, width = 0.01, height = 0.01) + theme_bw()
    grid.arrange(p1, p2)
    
    #Compute distance to optimal point : 1 recall and 1 PPV = full and pure cluster of targeted MoA
    res_cuts_perf$distToOptPoint = sqrt((1-res_cuts_perf$recall)**2 + (1-res_cuts_perf$ppv)**2)

    gc()
    return(res_cuts_perf)
}


