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

get_rect_positions <- function(ndosages, shift) {
  # we want to highlight the dosage with most interactions on the plot
  # Formula: 4 modes of action - the midpoints for the categories on the plot are 1, 2, 3, 4
  # The distance between groups of boxplots is 1/4
  # So the width of one boxplot is = (1-0.25) / (ndosages)
  # shift indicates how much the rectangle has to be shifted to indicate the desired dosage
  boxwidth <- 0.75/ndosages
  
  xmins <- rep(c(1 - boxwidth/2), 4) - ((ndosages %/% 2) * boxwidth)
  if ((ndosages %% 2) == 0) xmins <- xmins + boxwidth * 0.5 # with even number of dosages we have to nudge the rectangle to the right
  xmins <- xmins + shift * boxwidth
  xmins <- xmins + c(0, 1, 2, 3)
  
  xmaxs <- rep(c(1 + boxwidth/2), 4) - ((ndosages %/% 2) * boxwidth)
  if ((ndosages %% 2) == 0) xmaxs <- xmaxs + boxwidth * 0.5
  xmaxs <- xmaxs + shift * boxwidth
  xmaxs <- xmaxs + c(0, 1, 2, 3)
  
  ymins <- rep(-0.05, 4)
  ymaxs <- rep(1.05, 4)
  datafr <- data.frame(xmin = xmins, xmax = xmaxs,
    ymin = ymins, ymax = ymaxs)
  return(datafr)
}

ultimate_plot_new <- function(matrix_ext_row, filename, base_width = 100, base_height = 20) {
  # take a row from matrix_container_ext and plot probabilities for each drug for each repetition of the nested CV
  stopifnot(nrow(matrix_ext_row) == 1)
  # stopifnot("PredData" %in% names(matrix_ext_row))
  stopifnot("pred_data" %in% names(matrix_ext_row))
  
  predData <- matrix_ext_row$pred_data[[1]]
  # predData <- matrix_ext_row$PredData[[1]]
  
  plotlist <- map(unique(predData$process_broad), ~list())
  names(plotlist) <- unique(predData$process_broad)
  
  for (drug in unique(predData$drugname_typaslab)) {
    subfr <- predData[predData$drugname_typaslab == drug, ]
    my_moa <- unique(subfr$process_broad)
    stopifnot(length(my_moa) == 1)
    my_colours <- c("#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494")
    
    # get data frame for rectangle positions
    ndosages <- length(unique(subfr$conc))
    shift <- which(subfr$conc_mostias[seq_len(ndosages)]) - 1
    
    p <- ggplot(subfr, aes(x = moa_modelled, y = prob.moa, fill = factor(conc))) +
      # geom_rect(data = datafr, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = "grey50", alpha = 0.5, inherit.aes = FALSE) + # Error: Discrete value supplied to continuous scale - but putting at the end, geom_rect() works - what's going on???
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitterdodge(jitter.width = 0.2), shape = 1, size = 0.5) +
      geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_fill_manual(name = "Drug\nconc", values = my_colours) +
      labs(x = "MoA modelled", y = "Probability", title = paste0(drug, ": ", my_moa)) +
      geom_rect(data = get_rect_positions(ndosages = ndosages, shift = shift),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey50", alpha = 0.3, inherit.aes = FALSE)
    
    plotlist[[my_moa]] <- append(plotlist[[my_moa]], list(p))
  }
  
  # placeholder
  pblank <- ggplot(subfr, aes(x = moa_modelled, y = prob.moa, fill = factor(conc))) +
    geom_blank() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # fill up with empty plots for easier arrangement, hackhack
  for (elname in names(plotlist)) {
    el <- plotlist[[elname]]
    while(length(el) < max(lengths(plotlist))) {
      el <- append(el, list(pblank))
    }
    plotlist[[elname]] <- el
  }
  
  plotlist_flattened <- flatten(plotlist)
  ultimate_plot_new <- cowplot::plot_grid(plotlist = plotlist_flattened, ncol = max(lengths(plotlist)))
  
  cowplot::save_plot(filename, ultimate_plot_new, base_width = base_width, base_height = base_height, limitsize = FALSE)
}

plot_perf_from_container <- 
   function(containerObj, moa = c("cell_wall", "dna", "membrane_stress", "protein_synthesis"), 
            what = c("ROC", "prec-recall"), show_repeats = FALSE, row_var = NULL, col_var = NULL) {
      moa <- match.arg(moa, several.ok = TRUE)
      what <- match.arg(what)
      facet_flag <- TRUE 
      
      if (is.null(col_var) && (!is.null(row_var))) {
        stop("Please specify either row_var and col_var or only col_var but not only row_var")
      } else if (is.null(col_var) && is.null(row_var)) {
        if (nrow(containerObj) != 1) {stop("Without faceting nrow(containerObj) should be equal to 1")}
        facet_flag <- FALSE
      } else {
        facet_formula <- as.formula(paste(as.character(row_var), "~", as.character(col_var)))
        # check if number of subplots corresponds to number of rows:
        if (nrow(unique(model.frame(facet_formula, containerObj))) != nrow(containerObj)) {
          stop("Number of subplots won't match number of rows in containerObj")
        }
      }
      
      my_colours <- c("#e66101", "#fdb863", "#b2abd2", "#5e3c99")
      names(my_colours) <- c("cell_wall", "dna", "membrane_stress", "protein_synthesis")
      
      perf_measures <- 
        select(containerObj, perf_measures) %>%
        unnest() %>% 
        select(moa_modelled, auc, part_auc_01) %>%
        group_by(moa_modelled) %>%
        summarise(mean_auc = round(mean(auc), digits = 2), 
          mean_part_auc_01 = round(mean(part_auc_01), digits = 3))
      
      # make same order as in my_colours
      perf_measures <- perf_measures[match(perf_measures$moa_modelled, names(my_colours)), ]
      
      containerObj <- 
        select(containerObj, drug_dosages, feat_preselect, chemical_feats, fitted_model, 
          thresh_vs_perf) %>%
        unnest() %>%
        filter(moa_modelled %in% moa) %>%
        arrange(moa_modelled, cvrep, threshold)
      
      containerObj_averaged <-
         group_by(containerObj, moa_modelled, threshold, feat_preselect, fitted_model, drug_dosages, chemical_feats) %>%
         summarise(tpr_min = min(tpr), tpr_max = max(tpr), tpr_mean = mean(tpr),
           ppv_min = min(ppv), ppv_max = max(ppv), ppv_mean = mean(ppv), fpr_mean = mean(fpr)) %>%
         ungroup() %>%
         mutate(moa_modelled = factor(moa_modelled, levels = moa))
      
      p <- ggplot(containerObj_averaged, 
                  aes(x = if(what == "ROC") fpr_mean else tpr_mean, 
                      y = if(what == "ROC") tpr_mean else ppv_mean, 
                      group = moa_modelled, colour = moa_modelled, fill = moa_modelled))
      p <- p + geom_line(aes(colour = moa_modelled), size = 0.75) + 
      if (show_repeats) {
        p <- p + geom_path(data = containerObj, 
          aes(x = if(what == "ROC") fpr else tpr, y = if(what == "ROC") tpr else ppv, 
          group = interaction(cvrep, moa_modelled), colour = moa_modelled), alpha = 0.25, 
          inherit.aes = FALSE) 
        # side note: inherits.aes = FALSE is necessary here otherwise ggplot will complain about 
        # aesthetics from top level not present in containerObj (bug???)
      }
      if (facet_flag) p <- p + facet_grid(facet_formula)
      p <- p + 
        coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) + 
        labs(x = if(what == "ROC") "FPR (1-specificity)" else "TPR (recall)", 
          y = if(what == "ROC") "TPR (recall)" else "PPV (precision)") + 
        guides(fill = FALSE)
      
      if (!facet_flag & what == "ROC") {
        p <- p + 
        scale_colour_manual("Mode of action:", labels = paste0(names(my_colours), " \n(AUC: ", 
          perf_measures$mean_auc, "; AUC @10% FPR: ", perf_measures$mean_part_auc_01, ")\n"), 
          values = my_colours) + 
        theme(legend.position = c(0.7, 0.2))
      } else {
        p <- p + 
        scale_colour_manual("Mode of action:", labels = paste0(names(my_colours)), 
          values = my_colours) + 
        #theme(legend.position = c(0.7, 0.2))
        theme(legend.position = "bottom")
      }
      
      #return(list(plot = p, containerObj_averaged = containerObj_averaged, 
      #  containerObj = containerObj))
      return(p)
}

plot_prediction_probabilities <- 
  function(pred_data, thresh_vs_perf_data, moa, fileprefix, fpr_cutoff = 0.1) {
    # take a pred_data object (= from matrix_container, the PredData column) 
    # and plot prediction probabilities as a boxplot 
    # also add a thresh_vs_perf_data object to indicate the fpr cutoff 
    pred_data <- 
      filter(pred_data, moa_modelled == moa) %>% 
      select(prob.moa, drugname_typaslab, truth, conc)
    
    thresh_vs_perf_data <- 
      filter(thresh_vs_perf_data, moa_modelled == moa) %>%
      group_by(threshold) %>%
      summarise(mean_fpr = mean(fpr), mean_tpr = mean(tpr))
    
    # get threshold which is closest to an fpr of 0.1:
    thresh_vs_perf_data$diff <- abs(fpr_cutoff - thresh_vs_perf_data$mean_fpr)
    (thresh <- thresh_vs_perf_data$threshold[which.min(thresh_vs_perf_data$diff)])
    
    pred_data$drugndosg <- paste(pred_data$drugname_typaslab, pred_data$conc, sep = "_")
    pred_data$drugndosg <- fct_reorder(pred_data$drugndosg, .x = pred_data$prob.moa)
    
    selector <- 
      group_by(pred_data, drugname_typaslab, conc, drugndosg) %>%
      summarise(median_prob.moa = median(prob.moa)) %>%
      ungroup() %>%
      group_by(drugname_typaslab) %>%
      top_n(1, median_prob.moa) %>%
      pull(drugndosg)
    
    p <- ggplot(pred_data[pred_data$drugndosg %in% selector, ], aes(x = drugndosg, y = prob.moa)) + 
      geom_boxplot(aes(fill = truth), outlier.shape = 1, outlier.size = 1) + 
      geom_hline(yintercept = thresh, colour = "red") + 
      geom_hline(yintercept = 0.5, linetype = "dotted") + 
      coord_flip(ylim = c(0, 1)) + 
      theme_bw()
    ggsave(filename = paste0("./plots/", fileprefix, "_onedosg.pdf"), plot = p, width = 10, height = 15)
    
    # showing everything
    p <- ggplot(pred_data, aes(x = drugndosg, y = prob.moa)) + 
      geom_boxplot(aes(fill = truth), outlier.shape = 1, outlier.size = 1) + 
      geom_hline(yintercept = thresh, colour = "red") + 
      geom_hline(yintercept = 0.5, linetype = "dotted") + 
      coord_flip(ylim = c(0, 1)) + 
      theme_bw()
    ggsave(filename = paste0("./plots/", fileprefix, ".pdf"), plot = p, width = 10, height = 40)
  }


plot_highest_probs <- function(pred_data_obj, which_moa = c("cell_wall", "dna", "protein_synthesis", 
  "membrane_stress"), dosg_to_plot = c("highest_prob", "all")) 
{
  # takes pred_data table from matrix_container
  # have to specify the mode of action of interest
  # dosg_to_plot: either only plot the dosage of one drug that gets the highest probability 
  # or all dosages
  which_moa <- match.arg(which_moa)
  dosg_to_plot <- match.arg(dosg_to_plot)
  
  # pred_data_obj <- chosen_res$pred_data[[1]]
  pred_data_obj <- pred_data_obj[pred_data_obj$moa_modelled == which_moa, ]
  pred_data_obj$drug_conc <- paste(pred_data_obj$drugname_typaslab, pred_data_obj$conc, sep = "_")
  pred_data_obj$process_broad <- ifelse(pred_data_obj$process_broad == which_moa, 
    which_moa, paste0("not_", which_moa))
  stopifnot(nrow(pred_data_obj) > 0)
  
  winners <- 
    pred_data_obj %>%
    # first calculate median probability (across cv repeats) for each drug-concentration combination
    group_by(drugname_typaslab, drug_conc) %>%
    summarise(median_prob.moa = median(prob.moa)) %>%
    # then group them by drug and rank the median probabilities
    group_by(drugname_typaslab) %>%
    mutate(prob_rank = rank(median_prob.moa, ties.method = "first")) %>%
    # choose the drug-concentration combination with the highest rank
    arrange(desc(prob_rank), .by_group = TRUE) %>%
    top_n(1) %>% 
    arrange(desc(median_prob.moa))
  
  # and now keep in the original data (which isn't aggregated) the winning drug-concentration 
  # combination and define a factor so that all the data will be plotted in the correct order
  pred_data_obj_winners <- filter(pred_data_obj, drug_conc %in% winners$drug_conc)
  pred_data_obj_winners$drug_conc <- factor(pred_data_obj_winners$drug_conc, 
    levels = rev(winners$drug_conc))
  
  # ranking if all dosages are kept
  pred_data_obj$drug_conc <- forcats::fct_reorder(pred_data_obj$drug_conc, pred_data_obj$prob.moa)
  
  plot_func <- function(data) {
    ggplot(data, aes(x = drug_conc, y = prob.moa, fill = process_broad)) + 
      geom_boxplot(outlier.shape = 1, outlier.size = 0.5) + 
      coord_flip() + 
      labs(x = "Concentrations across all modes of action", y = "Drug-concentration combination", 
        title = paste0("Prediction probabilities across all CV repeats for ", which_moa))
  }
  
  if (dosg_to_plot == "highest_prob") {
    p <- plot_func(pred_data_obj_winners)
  } else {
    p <- plot_func(pred_data_obj)
  }
  
  print(p)
}
