get_outlier_boundaries <- function(data_vec) {
   # takes a numeric vector and returns the outlier boundaries according to boxplot.stats
   stopifnot(is.numeric(data_vec))
   outliers <- sort(boxplot.stats(data_vec)$out)
   lower_bound <- NA
   upper_bound <- NA
   if (any(outliers < 0)) {
      lower_bound <- max(outliers[outliers < 0])
   }
   if (any(outliers > 0)) {
      upper_bound <- min(outliers[outliers > 0])
   }
   
   data_vec_noout <- 
      if (is.na(lower_bound) & is.na(upper_bound)) { # no outliers = no boundaries
         data_vec
      } else if (!(is.na(lower_bound) | is.na(upper_bound))) { # i.e. we have 2 boundaries
         data_vec[data_vec > lower_bound & data_vec < upper_bound]
      } else if (is.na(lower_bound)) {
         data_vec[data_vec < upper_bound] 
      } else {
         data_vec[data_vec > lower_bound]
      }
   
   return(list(lower_bound = lower_bound, 
               upper_bound = upper_bound, 
               input_data = data_vec, 
               input_data_noout = data_vec_noout))
}
