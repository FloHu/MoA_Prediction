feature_distrib = function(m = the_matrix_top10pct, feature ){
    boxplot(m[, feature] ~ m$process_broad, lwd = 1.5, main = feature, col = rainbow(4), outline = F)
    # doesn't work - why?
    # stripchart(m[, feature] ~ m$process_broad, col = "black", bg = "grey", pch = 21, method = "jitter", add= T, vertical = T, cex = 1.5)
}
