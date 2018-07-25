make_importance_dfr <- function(feat_4model_output, nfeats = 20, normalise = TRUE, cutoff = 0.05) {
   # takes feature importances produced by plot_feat_4model() and turns them into a data frame
   # for easy plotting + labelling
   # normalise: to normalise highest importance to 1
   dfr <-  data.frame(importance = rev(feat_4model_output$features))
   dfr$gene <- factor(names(rev(feat_4model_output$features)),
                      levels = names(rev(feat_4model_output$features)))
   if (normalise) {
      dfr$importance <- dfr$importance / max(dfr$importance)
   }
   dfr$label <- ifelse(dfr$importance > cutoff, as.character(dfr$gene), "")
   return(dfr)
}

plot_importance <- function(dfr, moa) {
   # input: a data frame as produced by make_importance_dfr() and one MoA
   # output: a plot with feature importances, sorted in decreasing order. Leftmost point = most
   # important feature, normalised to one
   ggplot(dfr, aes(x = gene, y = importance)) +
      geom_point(shape = 1) +
      geom_text_repel(aes(label = label), size = 3) +
      #geom_text(aes(label = gene), size = 2, check_overlap = TRUE, nudge_x = 5) +
      theme_classic() +
      theme(axis.text.x = element_blank(), text = element_text(size = 16)) +
      labs(title = paste("Feature importances for MoA", moa), x = "Gene", y = "Importance (mean decrease in node impurity)")
}

feature_distrib_ggplot <- function(dfr, feature, save = TRUE) {
   # input: dfr is our original matrix (e.g. the_matrix_allDrugs_top10pct),
   # feature is a feature of your choice (e.g. RECA)
   # output is a number of boxplots how the corresponding feature is distributed split by MoA
   # Leo: this is the functional equivalent to your feature_distrib() function which I couldn't
   # get to work for printing to multiple pages as quickly as writing this function
   subdfr <- dfr[, c("drugname_typaslab", "process_broad", feature)]
   subdfr$process_broad <- factor(subdfr$process_broad,
                                  levels = c("dna", "cell_wall", "membrane_stress",
                                             "protein_synthesis", "oxidative_stress", "pmf",
                                             "protein_qc"))
   p <- ggplot(subdfr, aes_string(x = "process_broad", y = feature)) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(shape = 1, position = position_jitter(height = 0, width = 0.2)) +
      theme_classic() +
      theme(text = element_text(size = 16),
            axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12)) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      labs(x = "Mode of Action", y = paste0("s-score"),
           title = paste0("S-scores of all drugs in ", feature, " deletion (Nichols data)"))

   if (save) {
      ggsave(filename = paste0("./plots/present_feat_sscore_dist_", feature, ".pdf"), plot = p)
   }
   invisible(p)
}
