get_aucs_per_repeat <- function(predDataTbl) {
   #### TO DO: extract concatenated results objects already at the beginning, then they can be coerced to ROCR objects to get partial AUCs
   # take threshVsPerfData (= argument predDataTbl) and calculate auc for each repeat and each moa
   aucs <- 
      predDataTbl %>%
      mutate(negative = paste0("not_", moa_modelled)) %>%
      group_by(cvrep, moa_modelled) %>% # => we don't get aucs per fold but for the concatenation of all folds per repeat per moa
      summarise(auc = measureAUC(probabilities = prob.moa, truth = truth, negative = negative, positive = moa_modelled)) %>% 
      ungroup()
   return(aucs)
}
