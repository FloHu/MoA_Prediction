
distrib_drug_prob_allDosage = function(res, drug = "A22", dt_matrix = the_matrix_allDrugs){
    
    all_idDrug =  which(dt_matrix$drugname_typaslab == drug)
    drugMoa = as.character(dt_matrix[all_idDrug[1], "process_broad"] )
    
    par(mfrow = c(1,1))
    nbConc = length(all_idDrug)
    
    if(nbConc != 1){
        if(nbConc > 4){
            par(mfrow = c(3,2))
        }else if(nbConc ==2){
            par(mfrow = c(1,2))
        }else{
            par(mfrow = c(2,2))
        }
    }
    
    prob_allMoa = list()
    all_thres = list()
    
    for (idDrug in all_idDrug) {
        
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
        
        conc = as.character(dt_matrix[idDrug, "conc"] )
        boxplot(prob_allMoa, lwd = 2, col = rainbow(length(prob_allMoa)), ylim = c(0,1), outline = F,
                main = paste0(drug, " - ", drugMoa, " [C] = ", conc), cex.main = 1.5, cex.axis = 1, las = 1,
                sub = paste0("Dataset : ", deparse(substitute(res))))
        
        draw_thres_seg(all_thres)
        
        box(col = color_box_drug(all_thres, prob_allMoa, drugMoa), lwd = 10)
        stripchart(prob_allMoa, vertical = T, add = T, method = "jitter", pch = 21, bg = "grey", cex = 1.5)
        
    }
    
}
