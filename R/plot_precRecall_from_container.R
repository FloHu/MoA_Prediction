

plot_precRecall_from_container = function(containerObj, moa = "all", by = NULL){
    
    set_MoA = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
    match.arg(arg = moa, choices = c(set_MoA, "all"))
    
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
                summarise(ppv = mean(ppv), tpr = mean(tpr))
            
            plot_title = paste(containerObj[containerRow, ]$hyperparam_grid_name,
                               containerObj[containerRow, ]$drug_dosages,
                               containerObj[containerRow, ]$feat_preselect,
                               containerObj[containerRow, ]$chemical_feats, sep="_")
            
            p = ggplot(plotData, aes(x = tpr, y = ppv, colour = moa_modelled)) +
                geom_path(size = 1.5) +
                scale_color_manual(values = rainbow(length(set_MoA))) +
                labs(x = "True positive rate (Recall)", y = "Positive predictive value (Precision)", colour = "Mode of action", title = plot_title) +
                theme_bw()
            
            plot_list[[plot_title]] = p
        }
        
        #Ordering grid display
        if(!is.null(by)){
            filter1 = by[1]
            nrow_grid = length(table(containerObj[[filter1]]))
            filter2 = by[2]
            col_order = order(containerObj[[filter1]])
            plot_list = plot_list[col_order]
            
            if(!is.na(filter2)){
                ncol_grid = length(table(containerObj[[filter2]]))
                grid.arrange(arrangeGrob(grobs = plot_list, nrow = nrow_grid , ncol = ncol_grid,
                                         top =  textGrob(paste0("Columns ordered by ", by[2]) , gp=gpar(cex=2)),
                                         left = textGrob(paste0("Rows ordered by ", by[1]) , gp=gpar(cex=2), rot = 90)
                ))
                
            }else{
                ncol_grid = length(plot_list) / nrow_grid
                grid.arrange(arrangeGrob(grobs = plot_list, nrow = nrow_grid , ncol = ncol_grid, 
                                         left = textGrob(paste0("Rows ordered by ", by[1]) , gp=gpar(cex=2), rot = 90)
                ))
            }

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
                summarise(ppv = mean(ppv), tpr = mean(tpr))
            
            plotData$run_name  = paste(containerObj[containerRow, ]$hyperparam_grid_name,
                                       containerObj[containerRow, ]$drug_dosages,
                                       containerObj[containerRow, ]$feat_preselect,
                                       containerObj[containerRow, ]$chemical_feats,sep="_")
            plotDataMoa = rbind(plotDataMoa, plotData)
        }
        ggplot(plotDataMoa, aes(x = tpr, y = ppv, colour = run_name)) +
            geom_path(size = 1.5) +
            scale_color_manual(values = rainbow(nrow(containerObj))) +
            labs(x = "True positive rate (Recall)", y = "Positive predictive value (Precision)", colour = paste0("Model used _ MoA : ", moa )) +
            theme_bw()
        
    }
    
    
}






