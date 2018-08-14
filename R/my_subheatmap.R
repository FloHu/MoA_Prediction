my_subheatmap <- function(drugfeat_m, drugs, genes, ncolors = 70, palette = "cividis", from = -4, to = 4) {
   # take a drug_feature matrix (all numeric, conditions as rownames, features as column names) 
   # and a vector of drugs (condition) and gene names to get subportions of a heatmap (correlation-based clustering)
   library(viridis)
   library(circlize)
   my_cols <- match.fun(substitute(palette))(ncolors)
   
   stopifnot(is.numeric(drugfeat_m))
   stopifnot(is.matrix(drugfeat_m))
   
   if (any(range(nchar(colnames(drugfeat_m))) > 4)) {
      stop("Error: it looks like gene names are in the columns, not the rows.")
   }
   
   # get rows and columns
   drugpattern = paste(drugs, collapse = "|")
   genepattern = paste(genes, collapse = "|")
   drugfeat_m <- drugfeat_m[grepl(pattern = drugpattern, x = row.names(drugfeat_m)), 
                            grepl(pattern = genepattern, x = colnames(drugfeat_m)), 
                            drop = FALSE]
   
   if (!all(dim(drugfeat_m) >= 2)) {
      stop("Filtered matrix must have at least 2 rows in each dimension")
   }
   
   zerovarcols <- which(apply(drugfeat_m, 2, sd, na.rm = TRUE) == 0)
   if (any(zerovarcols)) {
      warning("Zero variance columns detected >:-( - not supported by correlation-based clustering.\n
           Continuing but with respective columns removed (", colnames(drugfeat_m)[zerovarcols], ")")
      drugfeat_m <- drugfeat_m[, -zerovarcols]
   }
   
   zerovarrows <- which(apply(drugfeat_m, 1, sd, na.rm = TRUE) == 0)
   if (any(zerovarrows)) {
      warning("Zero variance columns detected >:-( - not supported by correlation-based clustering.\n
           Continuing but with respective columns removed (", colnames(drugfeat_m[zerovarcols]), ")")
      drugfeat_m <- drugfeat_m[-zerovarrows, ]
   }
   
   # our distance function
   distfun <- function(x, y) {
      sqrt(1 - abs(cor(x, y)))
   }
   
   # to label cells
   cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.1f", drugfeat_m[i, j]), x, y, gp = gpar(col = "black", fontsize = 5))
   }
   
   # return(list(drugfeat_m, drugpattern, genepattern, from, to, ncolors, my_cols))
   
   Heatmap(drugfeat_m, col = colorRamp2(seq(from = from, to = to, length.out = ncolors), my_cols), 
           row_names_side = "left", 
           column_names_side = "top", 
           clustering_distance_rows = distfun, 
           clustering_distance_columns = distfun, 
           cell_fun = cell_fun)
}


