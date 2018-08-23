get_prob.moa_sds <- function(predDataTbl) {
   # take prediction data ("predDataTbl") and get for each drug the sd of the probabilities for which it is predicted 
   # to be a particular moa, once if moa_modelled == moa_actual and once moa_modelled != moa_actual 
   avg_sd_by_moa <- 
      predDataTbl %>%
      mutate(moa_modelled_is_truth = (moa_modelled == process_broad)) %>% # to distinguish the two cases 
      group_by(drugname_typaslab, moa_modelled_is_truth) %>% # we want to see stability for each drug depending on whether moa_modelled_is_truth 
      summarise(process_broad = unique(process_broad), 
                prob.moa_sd = sd(prob.moa), 
                prob.moa_sd2 = mean(tapply(prob.moa, conc, sd))) %>% # this is probably a better way to calculate the stability
      filter(process_broad != "other") # let's not consider these cases for now
   return(avg_sd_by_moa)
}
