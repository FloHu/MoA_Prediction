# all these functions work similarly: pass a nested cv instance + the corresponding data set

cvtest_check_n_obs <- function(nest_cvinst, data) {
   # checks if a cv instance has the same length as there are observations in the data
   # after that it checks if inner resampling instance has the same size as the size of the 
   # training set of the corresponding outer resampling instance
   required_length <- nrow(data)
   is_ok <- TRUE
   if (nest_cvinst$outer$size != required_length) {
      stop("Mismatch between cv instance size and number of data points")
      is_ok <- FALSE
   }
   
   if (!all(lengths(nest_cvinst$outer$train.inds)
            ==
            map_dbl(nest_cvinst$inner, "size"))) {
      stop("Mismatch between inner cv size and size of outer training set")
      is_ok <- FALSE
   }
   
   return(is_ok)
   # return(all(required_length == map_chr(cvinst, c("outer", "size"))))
}

cvtest_check_indices <- function(nest_cvinst, data) {
   ## make sure that no values occur twice in all the training/test splits
   # for outer instance:
   required_length <- nrow(data)
   is_ok <- TRUE
   
   # for each outer cv instance, make sure that no values occur twice
   inds <- map2(nest_cvinst$NCV_1$outer$train.inds, nest_cvinst$NCV_1$outer$test.inds, ~sort(c(.x, .y)))
   all(map_lgl(inds, function(.x) {
      if (!all(rle(.x)$lengths == 1) && all(rle(.x)$values) <= required_length) {
         stop("Error in train/test set splits (outer CV)")
         is_ok <- FALSE
      }
   }))
   
   # same for each inner cv instance
   inds <- map(nest_cvinst$inner, function(.x) {
      map2(.x[["train.inds"]], .x[["test.inds"]], ~sort(c(.x, .y)))
   })
   inds <- flatten(inds)
   if (!all(flatten_int(map(inds, function(.x) rle(.x)$lengths == 1)))) {
      stop("Error in train/test set splits of one of the inner CVs")
      is_ok <- FALSE
   }
   return(is_ok)
}

cvtest_check_blocking <- function(nest_cvinst, data) {
   # this function tests if blocking is respected
   
   # this is a list of indices that shouldn't be torn apart
   the_indivisible <- list()
   for (drugname in unique(data$drugname_typaslab)) {
      the_indivisible[[drugname]] <- which(drugname == data$drugname_typaslab)
   }
   
   check_exclusivity <- function(vec, train_test_pair) {
      # returns TRUE only if the elements of vec are found in either the first or the second element 
      # or in none of the elements of train_test pair, and FALSE if in both
      stopifnot(length(train_test_pair) == 2 || is.list(train_test_pair) || all(sapply(train_test_pair, is.numeric)))
      stopifnot(is.atomic(vec) || (!is.numeric(vec)))
      # think De Morgan!
      return(!all(vec %in% train_test_pair[[1]]) || !all(vec %in% train_test_pair[[2]]))
   }
   
   check_blocking <- function(cvinst, indexlist) {
      # function that "checks exclusivity" for a given CV instance
      # takes a list of indices and checks that in cvinst they are either in the training or the 
      # test set but not both
      is_ok <- TRUE
      train_test_pairs <- transpose(list(cvinst$train.inds, cvinst$test.inds))
      for (i in seq_along(indexlist)) {
         if (! all(map_lgl(train_test_pairs, check_exclusivity, vec = indexlist[[i]]))) {
            warning("An element was found in both the training and the test set (", names(indexlist)[i], ")!\n")
            is_ok <- FALSE
            break
         }
      }
      invisible(is_ok)
   }
   
   check_blocking(nest_cvinst$outer, indexlist = the_indivisible)
   
   # and now the same for the inner instances
   # first we need to translate the inner indices to the correct outer indices 
   # (inner indices correspond to indices of the vector "training.inds" of the corresponding outer 
   # training set)
   # within each inner instance: for each train/test pair needs to be matched up with the 
   # right outer training set indices
   
   # ensure we have the same amount of inner cv instances as training sets in the outer instance
   stopifnot(length(nest_cvinst$outer$train.inds) == length(nest_cvinst$inner))
   matching_numbers <- seq_len(length(nest_cvinst$outer$train.inds))
   
   for (m in matching_numbers) {
      outer_indices <- nest_cvinst$outer$train.inds[[m]]
      inner_inst <- nest_cvinst$inner[[m]]
      # map indices back
      inner_inst$train.inds <- map(inner_inst$train.inds, ~outer_indices[.x])
      inner_inst$test.inds <- map(inner_inst$test.inds, ~outer_indices[.x])
      
      # and now replace with 'correct' cv instance
      nest_cvinst$inner[[m]] <- inner_inst
   }
   
   # so now that the inner cv instances have been 'repaired' we can check_blocking() for each inner 
   # instance
   
   map_lgl(nest_cvinst$inner, check_blocking, indexlist = the_indivisible)
}

# this function is designed for visual inspection
cvtest_report_stratification <- function(nest_cvinst, data, return_as_table = TRUE) {
   # take a nested cv instance and corresponding data, then report how the different classes where 
   # split into training and test set
   # return_as_table: switch, if TRUE get back a table, otherwise the underlying data frame
   report_tables <- list()
   
   data$process_broad <- ifelse(data$process_broad %in% c("cell_wall", "dna", "membrane_stress", 
                                                          "protein_synthesis"), 
                                data$process_broad, 
                                "other")
   
   make_table <- function(cvinst, return_as_table) {
      train_test_pairs <- transpose(list(cvinst$train.inds, cvinst$test.inds), .names = letters[1:8])
      train_test_pairs <- map(train_test_pairs, function(.x) {names(.x) <- c('train', 'test'); return(.x)})
      train_test_pairs <- map(train_test_pairs, ~ unnest(enframe(.x)))
      
      train_test_pairs <- map2_dfr(train_test_pairs, names(train_test_pairs), function(.x, .y) {
         .x$split <- .y
         return(.x)
      })
      
      train_test_pairs$drugname <- data$drugname_typaslab[train_test_pairs$value]
      train_test_pairs$value <- data$process_broad[train_test_pairs$value]
      train_test_pairs <- select(train_test_pairs, value, name, split, drugname)
      names(train_test_pairs) <- c("MoA", "test.or.train", "split", "drugname")
      if (return_as_table) {
         train_test_pairs <- ftable(train_test_pairs[, c("MoA", "test.or.train", "split")])
      }
      return(train_test_pairs)
   }
   
   # translate inner indices to outer ones (see also 'cvtest_check_blocking')
   stopifnot(length(nest_cvinst$outer$train.inds) == length(nest_cvinst$inner))
   matching_numbers <- seq_len(length(nest_cvinst$outer$train.inds))
   
   for (m in matching_numbers) {
      outer_indices <- nest_cvinst$outer$train.inds[[m]]
      inner_inst <- nest_cvinst$inner[[m]]
      # map indices back
      inner_inst$train.inds <- map(inner_inst$train.inds, ~outer_indices[.x])
      inner_inst$test.inds <- map(inner_inst$test.inds, ~outer_indices[.x])
      
      # and now replace with 'correct' cv instance
      nest_cvinst$inner[[m]] <- inner_inst
   }
   
   report_tables[["outer"]] <- make_table(nest_cvinst$outer, return_as_table = return_as_table)
   report_tables[["inner"]] <- map(nest_cvinst$inner, make_table, return_as_table = return_as_table)
   return(report_tables)
}

cvtest_all_drugs_in_test <- function(nest_cvinst, data) {
   le_report <- cvtest_report_stratification(nest_cvinst, data, FALSE)
   
   check_table <- function(table) {
      all(unique(table$drugname) %in% table$drugname[table$test.or.train == "test"])
   }
   
   checks <- logical(length = 9)
   checks[1] <- check_table(le_report$outer)
   checks[2:9] <- map_lgl(le_report$inner, check_table)
   all(checks)
}

# master function
run_cv_tests <- function(rep_nest_cvinst, data) {
   cv_tests_report <- list()
   
   cv_tests_report$cvtest_check_n_obs <- 
      map_lgl(rep_nest_cvinst, cvtest_check_n_obs, data = data)
   
   cv_tests_report$cvtest_check_indices <- 
      map_lgl(rep_nest_cvinst, cvtest_check_indices, data = data)
   
   cv_tests_report$cvtest_check_blocking <- 
      map(rep_nest_cvinst, cvtest_check_blocking, data = data)
   
   cv_tests_report$cvtest_all_drugs_in_test <- 
      map_lgl(rep_nest_cvinst, cvtest_all_drugs_in_test, data = data)
   
   if (!all(unlist(cv_tests_report))) {
      cat("Potential problem with CV instance detected.")
      return(cv_tests_report)
   } else {
      invisible(TRUE)
   }
}

