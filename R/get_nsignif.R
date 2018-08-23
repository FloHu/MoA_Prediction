get_nsignif <- function(drugfeatmat, predDataTbl, the_matrix) {
   # function to add some more information to PredData: number of significant interactions, 
   # ranks of concentrations etc. 
   # drugfeatmat: from matrix_container_ext: the drug_feature matrix (used for modelling)
   # predDataTbl: from matrix_container_ext: the prediction table
   # the_matrix: from first notebook: the_matrix after the do(select_mutant(.)) step - contains 
   # sscores and qvalues for each drug-gene-concentration combination
   drugfeatmat_long <- 
      select(drugfeatmat, -process_broad) %>% 
      gather(3:ncol(.), key = "gene_synonym", value = "sscore")
   
   # by doing a join we can do a lookup which drug-concentration-gene combination is significant
   tab_nsignif <- left_join(drugfeatmat_long, the_matrix) %>%
      group_by(drugname_typaslab, conc) %>%
      summarise(n_signif = sum(significant, na.rm = TRUE))
   
   # and now join this back into the predDataTbl object - which is what we need for plotting
   # exclude drugs of category "other"
   tab_nsignif <- left_join(predDataTbl, tab_nsignif)
   tab_nsignif <- tab_nsignif[tab_nsignif$process_broad %in% 
                                 c("cell_wall", "dna", "membrane_stress", "protein_synthesis"), ]
   tab_nsignif$moa_modelled_is_truth <- tab_nsignif$moa_modelled == tab_nsignif$process_broad
   
   # in this case we're not interested in the probabilities for each repeat of the nested CV, only 
   # in the median
   tab_nsignif <- tab_nsignif %>%
      group_by(drugname_typaslab, conc, moa_modelled, moa_modelled_is_truth, n_signif, conc_mostias) %>%
      summarise(median_prob.moa = median(prob.moa))
   
   # get ranks of concentrations, of number of significant interactions, ...
   # but beware: conc with "highest number of significant interactions" was based on all features
   tab_nsignif <- 
      tab_nsignif %>%
      group_by(drugname_typaslab, moa_modelled) %>%
      mutate(conc_rank = rank(conc), 
             n_signif_rank = order(order(n_signif, conc)), # enables custom tie breaking as order(order(x)) == rank(x)
             n_conc = length(conc)) %>%
      ungroup(tmp)
   
   return(tab_nsignif)
}
