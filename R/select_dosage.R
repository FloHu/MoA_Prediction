select_dosage_most_ias <- function(dfr) {
   subfr <- select_if(ungroup(dfr), is.logical) # need to ungroup!
   sum_signif <- apply(subfr, 1, sum)
   # first occurrence wins, meaning the highest dosage will have the maximum rank in the case of ties:
   # (as dosages were sorted previously)
   ranks <- rank(sum_signif, ties.method = "first")
   ranks <- ranks == max(ranks)
   dfr$selector <- ranks
   return(dfr)
}
