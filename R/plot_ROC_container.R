


computeAUC = function(a){
    area = 0
    for(i in 1:(nrow(a)-1)){
        area = area + (a[i, ]$fpr - a[i+1, ]$fpr) * (a[i+1, ]$tpr + ((a[i, ]$tpr - a[i+1, ]$tpr) / 2))
    }
    return(area)
}



plot_ROC_from_container = function(containerObj, moa = "all", by = NULL){
    
    set_MoA = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    if(!moa %in% c(set_MoA, "all")){
        print("Invalid Mode of action")
        return(-1)
    }
   
    plot_list = list()
    
    library(grid)
    library(gridExtra)
    
    #================================= ALL MOAs =================================
    
    if(moa == "all"){
        for(containerRow in 1:nrow(containerObj)){
            
            perfData = containerObj[containerRow, ]$ThreshVsPerfData[[1]]
            predData = containerObj[containerRow, ]$PredData[[1]]
            
            perfData$threshold <- factor(perfData$threshold)
            # so it is possible to just average instead of concatenating folds
            plotData <- 
                perfData %>%
                group_by(threshold, moa_modelled) %>%
                summarise(fpr = mean(fpr), tpr = mean(tpr))
            
            #Adding AUC value correpsonding to the curves
            AUCdata <- list()
            for(m in set_MoA){
                AUCdata[[m]] = computeAUC(plotData %>% filter(moa_modelled == m))
            }
            
            AUCdata = round(unlist(AUCdata), digits = 3)
            plotData$moa_modelled = paste0(plotData$moa_modelled, " AUC : ",AUCdata[plotData$moa_modelled])
            
            plot_title = paste(containerObj[containerRow, ]$hyperparam_grid_name,
                               containerObj[containerRow, ]$drug_dosages,
                               containerObj[containerRow, ]$feat_preselect,
                               containerObj[containerRow, ]$chemical_feats, sep="_")
            
            p = ggplot(plotData, aes(x = fpr, y = tpr, colour = moa_modelled)) +
                geom_path(size = 1.5) +
                scale_color_manual(values = rainbow(length(set_MoA))) +
                geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5) +
                labs(x = "False positive rate", y = "True positive rate", colour = "Mode of action", title = plot_title) +
                theme_bw()
            
            plot_list[[plot_title]] = p
        }
        
        #Ordering grid display
        if(!is.null(by)){
            filter1 = by[1]
            nrow_grid = length(table(containerObj[[filter1]]))
            filter2 = by[2]
            if(!is.na(filter2)){
                ncol_grid = length(table(containerObj[[filter2]]))
            }else{
                ncol_grid = length(plot_list) / nrow_grid
            }
            col_order = order(containerObj[[filter1]])
            plot_list = plot_list[col_order]
            
            grid.arrange(arrangeGrob(grobs = plot_list, nrow = nrow_grid , ncol = ncol_grid))
        }else{
            grid.arrange(arrangeGrob(grobs = plot_list))
        }
        
    #================================= 1 MOA =================================    
        
    }else{
       
        plotDataMoa = NULL
        for(containerRow in 1:nrow(containerObj)){
            
            perfData = containerObj[containerRow, ]$ThreshVsPerfData[[1]]
            predData = containerObj[containerRow, ]$PredData[[1]]
            
            perfData = perfData %>% filter(moa_modelled == moa)
            
            perfData$threshold <- factor(perfData$threshold)
            # so it is possible to just average instead of concatenating folds
            plotData <- 
                perfData %>%
                group_by(threshold) %>%
                summarise(fpr = mean(fpr), tpr = mean(tpr))
            
            AUCdata = round(computeAUC(plotData), digits = 3)
            
            plotData$run_name  = paste(containerObj[containerRow, ]$hyperparam_grid_name,
                               containerObj[containerRow, ]$drug_dosages,
                               containerObj[containerRow, ]$feat_preselect,
                               containerObj[containerRow, ]$chemical_feats,
                               "AUC", AUCdata, sep="_")
            plotDataMoa = rbind(plotDataMoa, plotData)
        }
        ggplot(plotDataMoa, aes(x = fpr, y = tpr, colour = run_name)) +
            geom_path(size = 1.5) +
            scale_color_manual(values = rainbow(nrow(containerObj))) +
            geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", alpha = 0.5) +
            labs(x = "False positive rate", y = "True positive rate", colour = paste0("Model used _ MoA : ", moa )) +
            theme_bw()
        
    }

    
}






