# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.

# got it from https://gist.github.com/stevenworthington/3178163

ipak <- function(pkg){
  suppressMessages(suppressWarnings({
    pkg = as.character(substitute(pkg))
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    invisible(sapply(pkg, require, character.only = TRUE))
  }))
}
