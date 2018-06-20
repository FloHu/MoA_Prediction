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
  stopifnot("PredData" %in% names(matrix_ext_row))

  predData <- matrix_ext_row$PredData[[1]]

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
  ultimate_plot_new <- plot_grid(plotlist = plotlist_flattened, ncol = max(lengths(plotlist)))

  save_plot(filename, ultimate_plot_new, base_width = base_width, base_height = base_height, limitsize = FALSE)
}
