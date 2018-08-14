plot_pairplots <- function(genes, data = all_by_all) {
   # a very unflexible function that takes the all_by_all data frame and a vector with gene names 
   # to produce all possible combinations
   genes <- enquo(genes)
   
   filter(data, gene_synonym.x %in% !!genes, gene_synonym.y %in% !!genes) %>%
      ggplot(aes(x = s_score.x, y = s_score.y)) + 
      geom_point(aes(colour = process_broad), alpha = 0.5) + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      facet_wrap(gene_synonym.x ~ gene_synonym.y, strip.position = "bottom", scales = "free") + 
      theme_bw() + 
      labs(title = "Top strip = x-axis, bottom strip = y-axis")
}
