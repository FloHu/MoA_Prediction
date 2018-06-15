

draw_thres_seg = function(all_thres){
    segments(x0 = 0, x1 = 1.5, y0 = all_thres[[1]], lwd = 5, lty = "dotted")
    segments(x0 = 1.5, y0 = all_thres[[1]], y1 = all_thres[[2]],  lwd = 5, lty = "dotted")
    segments(x0 = 1.5, x1 = 2.5, y0 = all_thres[[2]], lwd = 5, lty = "dotted")
    segments(x0 = 2.5, y0 = all_thres[[2]], y1 = all_thres[[3]],  lwd = 5, lty = "dotted")
    segments(x0 = 2.5, x1 = 3.5, y0 = all_thres[[3]], lwd = 5, lty = "dotted")
    segments(x0 = 3.5, y0 = all_thres[[3]], y1 = all_thres[[4]],  lwd = 5, lty = "dotted")
    segments(x0 = 3.5, x1 = 5.5, y0 = all_thres[[4]], lwd = 5, lty = "dotted")
}


color_box_drug = function(all_thres, probMoa, drugMoa){

    med = lapply(probMoa, median)
    res = c()
    for (m in c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
        res = c(res, med[[paste0("prob_",m)]] > all_thres[[paste0("thres_",m)]])
    }
    names(res) = c("cell_wall", "dna", "membrane_stress", "protein_synthesis")

    if(sum(res) == 0){
        if(!drugMoa %in% c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
            return("green")
        }else{
            return("red")
        }
    }else if(sum(res) == 1){
        if(drugMoa %in% c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
            if(res[[drugMoa]]){
                return("green")
            }else{
                return("red")
            }
        }else{
            return("red")
        }
    }else{
        if(drugMoa %in% c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
            if(res[[drugMoa]]){
                return("orange")
            }else{
                return("red")
            }
        }else{
            return("red")
        }
    }
}


# MAIN FUNCTION
ultimate_plot = function(res, drugMat = the_matrix_allDrugs ){

    pdf(file = "~/Documents/ultimate_plot.pdf", height = 40, width = 200)
    #Defining layout
    
    eff = sort(table(drugMat$process_broad), decreasing = T)

    layout(mat = matrix(c(seq(1, eff[1]),
                            seq(from = eff[1]+1, length.out = eff[2]), rep(0, times = eff[1] - eff[2]),
                            seq(from = sum(eff[1:2])+1, length.out = eff[3]), rep(0, times = eff[1] - eff[3]),
                            seq(from = sum(eff[1:3])+1, length.out = eff[4]), rep(0, times = eff[1] - eff[4]),
                            seq(from = sum(eff[1:4])+1, to = sum(eff)), rep(0, 16)
                        ),
                    5, eff[1], byrow = TRUE
                    )
    )

    prob_allMoa = list()
    all_thres = list()
    
    # Generate Data to plot
    for (m in c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
        allData = c()
        for (r in 1:10){
            all_test_set = cat_outer_fold_pred(res = res, repetition = r, moa = m)
            allData = rbind(allData, all_test_set$data)
        }
        all_test_set$data = allData
        
        mccCurve = generateThreshVsPerfData(all_test_set, measures = mcc)
        #Best Threshold
        bestThres = mccCurve$data[which.max(mccCurve$data$mcc),"threshold"]
        all_thres[[paste0("thres_", m)]] = bestThres
        prob_allMoa[[paste0("prob_", m)]] = all_test_set
    }
    
    drugMat$process_broad = factor(drugMat$process_broad, levels = names(eff))
    
    #Filter data to plot for each drug and plot it
    for (d in drugMat[order(drugMat$process_broad), "drugname_typaslab"]) {
        idDrug =  which(drugMat$drugname_typaslab == d)
        drugMoa = as.character(drugMat[idDrug, "process_broad"])
       
        plotData_drug = list()
        for (m in c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
            plotData_drug[[paste0("prob_", m)]] = prob_allMoa[[paste0("prob_", m)]]$data %>% filter( id == idDrug ) %>% select(paste0("prob.", m)) %>% t(.)
            
        }
        boxplot(plotData_drug, lwd = 2, col = rainbow(length(plotData_drug)), ylim = c(0,1), outline = F,
                main = paste0(d, " - ", drugMoa), cex.main = 3, cex.axis = 2, las = 1)
        
        draw_thres_seg(all_thres)
        box(col = color_box_drug(all_thres, probMoa = plotData_drug, drugMoa), lwd = 10)
        stripchart(plotData_drug, vertical = T, add = T, method = "jitter", pch = 21, bg = "grey", cex = 1.5)
    }
    
    dev.off()
}


# ==============================================================================

distrib_drug_prob_small = function(res, drug = "A22", dt_matrix = the_matrix_allDrugs){

    idDrug =  which(dt_matrix$drugname_typaslab == drug)
    drugMoa = as.character(filter(dt_matrix, drugname_typaslab == drug) %>% select(process_broad))

    prob_allMoa = list()
    all_thres = list()

    for (m in c("cell_wall", "dna", "membrane_stress", "protein_synthesis")){
        allData = c()
        for (r in 1:10){
            all_test_set = cat_outer_fold_pred(res = res, repetition = r, moa = m)
            allData = rbind(allData, all_test_set$data)
        }
        all_test_set$data = allData

        mccCurve = generateThreshVsPerfData(all_test_set, measures = mcc)
        #Best Threshold
        bestThres = mccCurve$data[which.max(mccCurve$data$mcc),"threshold"]
        all_thres[[paste0("thres_", m)]] = bestThres
        prob_allMoa[[paste0("prob_", m)]] = filter(all_test_set$data, id == idDrug ) %>% select(paste0("prob.", m)) %>% t(.)
    }


    boxplot(prob_allMoa, lwd = 2, col = rainbow(length(prob_allMoa)), ylim = c(0,1), outline = F,
            main = paste0(drug, " - ", drugMoa), cex.main = 1.5, cex.axis = 1, las = 1,
            sub = paste0("Dataset : ", deparse(substitute(res))))

    draw_thres_seg(all_thres)

    box(col = color_box_drug(all_thres, prob_allMoa, drugMoa), lwd = 10)

    stripchart(prob_allMoa, vertical = T, add = T, method = "jitter", pch = 21, bg = "grey", cex = 1.5)

}
