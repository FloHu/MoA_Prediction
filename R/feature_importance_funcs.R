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

plot_heatmap <- function(m = the_matrix_allDrugs_top10pct, feats_to_keep = "all", filename = "heatmap.pdf") {
   # plots a heatmap for all drugs with only the features provided in the feats_to_keep argument
   cluster_matrix <- m
   row.names(cluster_matrix) <- cluster_matrix$drugname_typaslab
   cluster_matrix$drugname_typaslab <- NULL
   process_lut <- cluster_matrix$process_broad
   names(process_lut) <- row.names(cluster_matrix)
   cluster_matrix$process_broad <- NULL
   cluster_matrix <- as.matrix(cluster_matrix)

   # keep only genes in the columns because of measurement scale
   if (feats_to_keep != "all") {
      cluster_matrix <- cluster_matrix[, colnames(cluster_matrix) %in% feats_to_keep]
   }

   # get correlation-based distances between the drugs (--> transpose matrix)
   drug_drug_dist <- as.dist(1 - abs(cor(t(cluster_matrix))))
   gene_gene_dist <- as.dist(1 - abs(cor(cluster_matrix)))

   # define colors for the mode of action
   moa_to_colour <- c(cell_wall = "#a6cee3",
                      dna = "#1f78b4",
                      protein_synthesis = "#b2df8a",
                      membrane_stress = "#33a02c",
                      pmf = "#bababa", protein_qc = "#bababa", oxidative_stress = "#bababa"
                      )
   rowcols <- moa_to_colour[process_lut[rownames(cluster_matrix)]]

   pdf(paste0("./plots/", filename), width = 20, height = 20)
   heatmap.2(cluster_matrix,
             trace = "none",
             breaks = c(-10, -5, -3, -1, 0, 1, 3, 5, 10), # use this i/o zlim
             Rowv = as.dendrogram(hclust(drug_drug_dist)), # change order of rows
             Colv = as.dendrogram(hclust(gene_gene_dist)),
             RowSideColors = rowcols,
             margins = c(12, 9),
             dendrogram = "row") # only draw row dendrogram
   legend("top",
          legend = names(moa_to_colour),
          col = moa_to_colour,
          lty = 1,
          lwd = 10)
   dev.off()
}
