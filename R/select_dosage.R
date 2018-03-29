select_dosage <- function(dfr) {
   if ( all(dfr$significant) || all(!dfr$significant) ) {
      return(dfr[which.max(dfr$conc), ])
   } else {
      sel_conc <- max(dfr$conc[dfr$significant])
      return(dfr[dfr$conc == sel_conc, ])
   }
}
