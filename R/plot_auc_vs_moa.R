plot_auc_vs_moa <- function(dataset, xaxis, yaxis, colourvar, 
                            facet_expr = quote(facet_grid(fitted_model ~ .)), 
                            additional_info, showPlot = TRUE) {
   p <- 
      ggplot(dataset, aes_string(x = xaxis, y = yaxis, colour = colourvar)) + 
      geom_boxplot(position = position_dodge(), outlier.shape = NA) + 
      geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 1, alpha = 0.5) + 
      eval(facet_expr) + 
      labs(x = "Mode of action modelled", y = paste0(yaxis, "(one dot per repeat of nested CV)"), 
           title = paste0("Performances of MoA predictions\n", additional_info)) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 18))
   
   if (showPlot) print(p)
   invisible(p)
}
