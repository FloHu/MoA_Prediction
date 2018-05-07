plot_dend = function(drugs_fingerprint, ...){
    # purpose: pass matrices with chemical fingerprints (similarities between drug pairs)
    # and plot with hierarchical clustering to get similarity
    fp_matrix = apply(drugs_fingerprint, 1,  function(x) {
                unlist(strsplit(x["RDkit_fingerprint"], split = "" ))
        }
    )
    fp_matrix = t(fp_matrix)
    rownames(fp_matrix) = drugs_fingerprint$drugnames_typaslab

    color_moa = c(rainbow(4), rep("black", 3))
    names(color_moa) = c("dna", "cell_wall", "membrane_stress", "protein_synthesis", "protein_qc", "oxidative_stress", "pmf")

    dist_fp_matrix = dist(fp_matrix, method = "binary")

    hclust_fp = hclust(dist_fp_matrix)
    dend_fp = as.dendrogram(hclust_fp, hang = 0.1)

    labels_colors(dend_fp) <- color_moa[the_matrix_allDrugs[order.dendrogram(dend_fp), "process_broad"]]
    par(mar = c(3,2,2,8))
    dend_fp %>%
        set("labels_cex", 0.7) %>%
        set("branches_lwd", 1.5) %>%
        plot(horiz = T, cex = 0.5, ...)

    abline(v = 0.1, lty = "dotted", lwd = 2, col ="darkgrey")
}
