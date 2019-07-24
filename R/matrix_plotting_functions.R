plot_auc_vs_moa <- function(dataset, xaxis, yaxis, colourvar, groupvar,
  facet_expr = quote(facet_grid(fitted_model ~ .)),
  xlab, ylab, title, additional_info, showPlot = TRUE) {
  p <- ggplot(dataset, aes_string(x = xaxis, y = yaxis, colour = colourvar,
    group = groupvar)) +
    # geom_boxplot(position = position_dodge(), outlier.shape = NA) +
    # geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 1, alpha = 0.5) +
    stat_summary(geom = "line", fun.y = median) +
    stat_summary(geom = "point", fun.y = median) +
    stat_summary(geom = "errorbar", fun.y = median, width = 0.2, alpha = 0.5) +
    eval(facet_expr) +
    labs(x = xlab, y = ylab, title = title, subtitle = additional_info) +
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
    while (length(el) < max(lengths(plotlist))) {
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
            what = c("ROC", "prec-recall"), show_repeats = FALSE, row_var = NULL, col_var = NULL)
     {
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

      my_colours <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
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
         group_by(containerObj, moa_modelled, threshold, feat_preselect,
           fitted_model, drug_dosages, chemical_feats) %>%
         summarise(tpr_min = min(tpr), tpr_max = max(tpr), tpr_mean = mean(tpr),
           ppv_min = min(ppv), ppv_max = max(ppv), ppv_mean = mean(ppv),
           fpr_mean = mean(fpr)) %>%
         ungroup() %>%
         mutate(moa_modelled = factor(moa_modelled, levels = moa))

      p <- ggplot(containerObj_averaged,
                  aes(x = if (what == "ROC") fpr_mean else tpr_mean,
                      y = if (what == "ROC") tpr_mean else ppv_mean,
                      group = moa_modelled, colour = moa_modelled,
                    fill = moa_modelled))
      p <- p + geom_line(aes(colour = moa_modelled), size = 0.75) +
      if (show_repeats) {
        p <- p + geom_path(data = containerObj,
          aes(x = if (what == "ROC") fpr else tpr, y = if (what == "ROC") tpr else ppv,
          group = interaction(cvrep, moa_modelled), colour = moa_modelled), alpha = 0.25,
          inherit.aes = FALSE)
        # side note: inherits.aes = FALSE is necessary here otherwise ggplot will complain about
        # aesthetics from top level not present in containerObj (bug???)
      }
      if (facet_flag) p <- p + facet_grid(facet_formula)
      p <- p +
        coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
        labs(x = if (what == "ROC") "FPR (1-specificity)" else "TPR (recall)",
          y = if (what == "ROC") "TPR (recall)" else "PPV (precision)") +
        guides(fill = FALSE)

      if (!facet_flag & what == "ROC") {
        p <- p +
        scale_colour_manual("Mode of action:", labels = paste0(names(my_colours), " \nAUC: ",
          perf_measures$mean_auc, "\nAUC @10% FPR: ", perf_measures$mean_part_auc_01, "\n"),
          values = my_colours) #+
        #theme(legend.position = c(0.7, 0.3))
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

plot_mcl_probs_lines <- function(melted_pred_data, printplot = TRUE, 
  save = FALSE, file = NULL) {
  # needs a melted multiclass prediction (function melt_pred_data)
  p <- ggplot(melted_pred_data, aes(x = factor(conc), y = prob.med,
    colour = predicted_prob)) +
    geom_linerange(aes(ymin = prob.min, ymax = prob.max), size = 0.5) +
    geom_line(aes(group = predicted_prob), size = 0.25) + 
    geom_point(size = 0.5) + 
    facet_wrap( ~ truth + drugname_typaslab, scales = "free_x") +
    geom_hline(yintercept = c(0.25), linetype = c("dotted")) +
    scale_colour_manual("Predicted probability\nfor class",
      values = moa_cols2, labels = moa_repl2) +
    coord_cartesian(ylim = c(0, 1)) + 
    theme(text = element_text(size = 5), panel.grid = element_blank()) + 
    labs(x = "Concentration", y = "Probability")

  if (printplot) print(p)
  
  if (save) {
    if (is.null(file)) stop("Must provide a filename")
    suppressWarnings(
      ggsave(filename = file, plot = p, width = 180, height = 200, units = "mm")
    )
  }
  
  invisible(p)
}

plot_mcl_probs_heatmap <- function(melted_pred_data, mics, printplot = TRUE, 
  save = FALSE, file = NULL) {
  tmp <- suppressMessages(left_join(melted_pred_data, mics))
  tmp <- tmp %>%
    mutate(drug_conc = sprintf("%-19s %-5s %-5s", 
      Hmisc::capitalize(tolower(drugname_typaslab)), conc, mic_curated),
      prob.med.range = cut(prob.med, breaks = seq(from = 0, to = 1, by = 0.1),
        labels = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%",
          "60-70%", "70-80%", "80-90%", "90-100%")),
      drug_conc = factor(drug_conc,
        levels = unique(drug_conc[order(truth, drugname_typaslab, conc)])))

  tmp$truth <- fct_recode(tmp$truth, "Annotated MoA: Cell Wall" = "cell_wall",
    "Annotated MoA: DNA" = "dna", "Annotated MoA: Membrane Stress" = "membrane_stress",
    "Annotated MoA: Protein Synthesis" = "protein_synthesis")

  tmp <- group_by(tmp, drug_conc) %>%
    mutate(is_max = (prob.max) == max(prob.max))
  tmp$geompoint <- ifelse(tmp$is_max, tmp$drug_conc, NA)
  tmp$geompoint <- levels(tmp$drug_conc)[tmp$geompoint]

  # cols <- RColorBrewer::brewer.pal(9, "BuPu")
  # cols <- colorRampPalette(cols)(10)
  cols <- RColorBrewer::brewer.pal(min(9, nlevels(tmp$prob.med.range)), "BuPu")
  cols <- colorRampPalette(cols)(min(10, nlevels(tmp$prob.med.range)))
  
  tmp <- ungroup(tmp)
  tmp$drug_conc %<>% droplevels()
  tmp$drug_conc <- fct_relevel(tmp$drug_conc, rev(levels(tmp$drug_conc)))

  p <- ggplot(tmp, aes(x = drug_conc, y = predicted_prob)) +
    geom_tile(aes(fill = prob.med.range)) +
    geom_point(aes(x = geompoint), size = 0.3) +
    #coord_flip() +
    scale_fill_manual("Probability", values = cols) +
    scale_y_discrete("MoA predicted", labels = moa_repl2) +
    scale_x_discrete("") + 
    # scale_x_discrete("Drug + concentration", limits = levels(tmp$drug_conc)) + 
    paper_theme + 
    theme(axis.text.y = element_text(family = "Courier", size = 6), 
      legend.position = "bottom") + 
    coord_flip() + 
    facet_wrap( ~ truth, drop = TRUE, scales = "free")
  p
  
  if (printplot) suppressWarnings(print(p))

  if (save) {
    if (is.null(file)) stop("Must provide a filename")
    suppressWarnings(
      ggsave(filename = file, plot = p, width = 190, height = 250, units = "mm")
    )
  }
  
  invisible(p)
}

plot_mcl_probs_heatmap2 <- function(pred_data, 
  mics, 
  printplot = TRUE, 
  save = FALSE, 
  file = NULL, 
  order_by_conf = FALSE) {
  
  tmp <- suppressMessages(left_join(pred_data, mics))
  tmp$drug_conc = paste0(tmp$drugname_typaslab, "_", tmp$conc, " (", 
    tmp$mic_curated, ")")
  
  tmp <- modify_if(tmp, is.factor, as.character) %>%
    group_by(drugname_typaslab, conc, drug_conc) %>%
    summarise(prob.cell_wall = median(prob.cell_wall), 
      prob.dna = median(prob.dna), 
      prob.membrane_stress = median(prob.membrane_stress), 
      prob.protein_synthesis = median(prob.protein_synthesis), 
      truth = unique(truth)) %>%
    gather(prob.cell_wall:prob.protein_synthesis, key = "MoA", value = "prob") %>%
    arrange(truth, drugname_typaslab, conc) %>% 
    mutate(MoA = str_replace(string = MoA, pattern = "prob\\.", replacement = "")) %>%
    ungroup()
  
  # now derive a response and a delta-prob (as confidence score)
  tmp <- group_by(tmp, drugname_typaslab, conc) %>%
    mutate(prob.delta = max(prob) - max(prob[-which.max(prob)])) %>% 
    slice(which.max(prob)) %>%
    rename(response = MoA) %>%
    ungroup()
  
  tmp <- tmp %>%
    mutate(prob.delta.range = cut(prob.delta, breaks = seq(from = 0, to = 1, by = 0.1), 
      labels = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", 
        "60-70%", "70-80%", "80-90%", "90-100%")), 
      drug_conc = factor(drug_conc, 
        levels = unique(drug_conc[order(truth, drugname_typaslab, conc)])))
  
  tmp$truth <- fct_recode(tmp$truth, "Annotated MoA: Cell Wall" = "cell_wall",
    "Annotated MoA: DNA" = "dna", "Annotated MoA: Membrane Stress" = "membrane_stress",
    "Annotated MoA: Protein Synthesis" = "protein_synthesis")
  
  if (order_by_conf) {
    tmp$drug_conc <- fct_reorder(tmp$drug_conc, tmp$prob.delta)
  }
  
  p <- ggplot(tmp, aes(x = drug_conc, y = prob.delta)) +
    geom_bar(aes(fill = response), stat = "identity", width = 0.75) +
    scale_fill_manual("Response\n(MoA predicted)", values = moa_cols, labels = moa_repl) +
    scale_y_continuous("Confidence of prediction: delta of probabilities: highest vs. 2nd-highest") +
    scale_x_discrete("Drug + concentration") + 
    coord_flip() + 
    comparison_theme + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
      panel.grid = element_blank()) + 
    facet_wrap( ~ truth, scales = "free_y")
  
  if (printplot) suppressWarnings(print(p))
  
  if (save) {
    if (is.null(file)) stop("Must provide a filename")
    suppressWarnings(
      ggsave(filename = file, plot = p, width = 180, height = 200, units = "mm")
    )
  }
  
  invisible(p)
}

plot_perf <- function(mc_ext, what, save = FALSE, file = NULL) {
  # input: matrix container (mc_ext); what = which performance measure to plot
  # (e.g. mmce, kappa)
  # output: a plot with the different algorithms on the x-axis, chemical_feats
  # and drug_dosages displayed by facets

  p <-
  select(mc_ext, fitted_model, drug_dosages, chemical_feats, perf_measures) %>%
  unnest() %>%
  group_by(fitted_model, drug_dosages, chemical_feats, cvrep) %>%
  summarise(mmce_mean = mean(mmce), kappa_mean = mean(kappa)) %>%
  ggplot(aes_string(x = "fitted_model", y = what)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.25, height = 0),
      shape = 1) +
    # need to change labeller function in facet_grid to change faceting labels,
    # see also
    # https://stackoverflow.com/questions/3472980/how-to-change-facet-labels
    # and the ggplot documentation
    facet_grid(chemical_feats ~ drug_dosages) +
    labs(y = paste0("Performance (", what, ")")) +
    scale_x_discrete("Classifier type", labels = classifier_repl) +
    comparison_theme + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)

  if (save) {
    if (is.null(file)) stop("Must provide a filename")
    ggsave(filename = file, plot = p, width = 87, height = 80, units = "mm")
  }
}

plot_inner_vs_outer <- function(mc_ext, save = FALSE, file = NULL) {
  # are performances of inner cross-validations similar to performances of 
  # the test set? 
  # input: matrix container (mc_ext)
  # output: a plot comparing mmce and kappa between inner cross-validated 
  # performance measure and test set performance measure
  tmp <- mc_ext %>%
    select(fitted_model, drug_dosages, chemical_feats, perf_measures, 
      opt_pars) %>%
    unnest()
  stopifnot({
    all(tmp$cvrep == tmp$cvrep1)
    all(tmp$split == tmp$split1)
  })
  
  p1 <- ggplot(tmp, aes(x = mmce.cv.inner, y = mmce, colour = cvrep)) + 
    geom_point() + 
    coord_fixed() + 
    comparison_theme
  
  p2 <- ggplot(tmp, aes(x = "", y = mmce - mmce.cv.inner)) + 
    geom_boxplot(outlier.shape = NA, width = 0.5) + 
    ggbeeswarm::geom_beeswarm(cex = 3, aes(colour = cvrep)) + 
    comparison_theme
  
  p3 <- ggplot(tmp, aes(x = kappa.cv.inner, y = mmce, colour = cvrep)) + 
    geom_point() + 
    coord_fixed() + 
    comparison_theme
  
  p4 <- ggplot(tmp, aes(x = "", y = kappa - kappa.cv.inner)) + 
    geom_boxplot(outlier.shape = NA) + 
    ggbeeswarm::geom_beeswarm(cex = 3, aes(colour = cvrep)) + 
    comparison_theme
  
  # tutorial: 
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
  p <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, 
    nrow = 2)
  print(p)
  
  if (save) {
    if (is.null(file)) stop("Must provide a filename")
    ggsave(filename = file, plot = p, width = 180, height = 180, units = "mm")
  }
}

plot_heatmap <- function(dfm, feats, mics, moa, split = NULL, printplot = FALSE, save = FALSE, 
  file = NULL) {
  # dfm = drug-feature matrix
  # feats = features to keep for dfm
  # mics = info on mics
  # moa = table about MoAs
  m <- dfm
  m <- suppressMessages(
    left_join(m, moa[, c("drugname_typaslab", "process_subgroup")]) %>%
    left_join(mics[, c("drugname_typaslab", "mic_curated")])
  )
  m <- arrange(m, process_broad, process_subgroup, drugname_typaslab)
  # row_order <- m$drugname_typaslab
  
  dfr <- as.data.frame(m[, c("process_broad")])
  rownames(m) <- paste0(m$drugname_typaslab, "_", m$conc, " (", 
    m$mic_curated, ")")
  
  m <- select(m, -drugname_typaslab, -conc, -process_broad, -process_subgroup, 
    -mic_curated) %>%
    select(feats) %>% # cannot be mixed with previous select statement
    as.matrix()
  # m[1:5, 1:5]
  
  # Make annotation object
  # ComplexHeatmap seems to have problems with tibbles! 
  row_annot <- HeatmapAnnotation(df = dfr, 
    col = list(process_broad = moa_cols), which = "row", 
    name = "MoA", annotation_width = 3)
  # row_annot
  # draw(row_annot, 1:20)
  
  my_distfun <- function(x, y) {1 - abs(cor(x, y))}
  min_col <- plasma(2)[1]
  max_col <- plasma(2)[2]
  
  h <- Heatmap(matrix = m, 
    col = colorRamp2(breaks = c(-5, 5), colors = c(min_col, max_col)), 
    name = "S-score", 
    clustering_distance_rows = my_distfun, 
    clustering_distance_columns = my_distfun, 
    column_names_side = "top", 
    row_names_side = "right", 
    cluster_rows = TRUE, 
    row_dend_side = "right", 
    row_names_gp = gpar(fontsize = 4), 
    column_names_gp = gpar(fontsize = 5), 
    split = split, 
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", m[i, j]), x, y, 
        gp = gpar(col = "black", fontsize = 2))
    }) + row_annot
  
  if (printplot) print(h)
  
  if (save) {
    if (is.null(file)) stop("Must provide a filename")
    pdf(file = file, width = 7, height = 15)
    print(h)
    dev.off()
  }
  
  invisible(h)
}

plot_tsne <- function(tsne_mat, seeed = 5, dims = 3, perplexity = 10, max_iter = 6000, 
  save = FALSE, file = NULL) {
  tsne_mat$process_broad <- factor(tsne_mat$process_broad)
  set.seed(seeed)
  
  stopifnot(all(colnames(tsne_mat)[1:3] %in% c("drugname_typaslab", "conc", "process_broad")))
  
  tsne_res <- Rtsne::Rtsne(X = tsne_mat[, -c(1:3)], dims = dims, perplexity = perplexity, 
    max_iter = max_iter)
  
  plotData <- data.frame(tsne_res$Y, MoA = tsne_mat$process_broad, 
    drug = paste0(tsne_mat$drugname_typaslab, "_", tsne_mat$conc))
  colnames(plotData)[1:3] <- c("tSNE1", "tSNE2", "tSNE3")
  
  p <- plot_ly(plotData, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, color = ~MoA, 
    colors = moa_cols[levels(tsne_mat$process_broad)], 
    marker = list(size = 8, line = list(color = 'rgba(0, 0, 0, 1)', 
      width = 1.5))) %>%
    add_markers(text = ~ drug) %>%
    layout(title = "tSNE map", 
      scene = list(xaxis = list(title = 'tSNE1'), yaxis = list(title = 'tSNE2'), 
        zaxis = list(title = 'tSNE3')))
  print(p)
  
  if (save) {
    if (is.null(file)) stop()
    saveRDS(object = p, file = file)
  }
  
  plotlist <- list()
  p1 <- ggplot(plotData, aes(x = tSNE1, y = tSNE2, colour = MoA)) + 
    geom_point(alpha = 0.75) + 
    scale_colour_manual("MoA", values = moa_cols, labels = moa_repl) + 
    comparison_theme
  plotlist[["p1"]] <- p1
  p2 <- p1 + aes(y = tSNE3)
  plotlist[["p2"]] <- p2
  p3 <- p2 + aes(x = tSNE2)
  plotlist[["p3"]] <- p3
  
  if (save) {
    pdf(file = paste0("./plots/", basename(file), "_2D.pdf"), width = 3.42, height = 2.8)
    walk(plotlist, print)
    dev.off()
  }
}

plot_pca <- function(mat, save = FALSE, file = NULL) {
  # input: drug-feature matrix (wide format)
  
  # remove "metadata"
  m <- select(mat, -one_of(c("drugname_typaslab", "conc", "process_broad")))
  stopifnot(every(m, is.numeric))
  
  rownames(m) <- paste0(mat$drugname_typaslab, "_", mat$conc)
  pr.out <- prcomp(m, scale = TRUE)
  # for screeplot
  pve <- pr.out$sdev ^ 2
  pve <- pve / sum(pve)
  names(pve) <- paste0("PC", seq_along(pve))
  pve <- enframe(pve)
  pve$name <- factor(pve$name, levels = pve$name)
  # 4 pages of plots: first page: screeplot, other pages: pairwise combinations 
  # of principal components
  screepl <- ggplot(pve[1:30, ], aes(x = name, y = value)) + 
    geom_path(aes(group = 1)) + 
    geom_point(shape = 1) + 
    comparison_theme + 
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(title = "Scree plot of first 30 PCs", x = "Principal component", 
      y = "Fraction of variance explained")
  
  prcomps <- as.data.frame(pr.out$x)
  prcomps$drug_conc <- row.names(prcomps)
  prcomps <- separate(prcomps, col = drug_conc, 
    into = c("drugname_typaslab", "conc"), sep = "_")
  prcomps <- left_join(prcomps, unique(mat[, c("drugname_typaslab", "process_broad")]))
  prcomps <- select(prcomps, drugname_typaslab, conc, process_broad, everything())
  
  my_legend <- scale_colour_manual("Target process", labels = moa_repl, 
    values = moa_cols)
  
  # there are some interesting dots in the PC1-PC2 component that I would like to 
  # label
  
  p1 <- ggplot(prcomps, aes(x = PC1, y = PC2)) + 
    geom_point(aes(colour = process_broad), alpha = 0.75) + 
    theme(panel.grid = element_blank(), legend.position = "top") + 
    comparison_theme + 
    my_legend
  p2 <- p1 + aes(x = PC1, y = PC3) + theme(legend.position = "None")
  p3 <- p2 + aes(x = PC2, y = PC3)
  
  p <- grid.arrange(screepl, p1, p2, p3, nrow = 2)
  print(p)
  
  if (save) {
    if (is.null(file)) stop("Must provide a filename")
    ggsave(filename = file, plot = p, width = 180, height = 160, units = "mm")
  }
  
  invisible(p)
}



