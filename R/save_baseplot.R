save_baseplot <- function(filename = "./plots/save_baseplot.pdf", plot, width = 7, height = 7) {
  # takes a filename, derives the graphics device from the extension and saves the provided recordedplot (function recordPlot()) to a file
  library(stringr)
  device <- match.fun(str_extract(filename, pattern = "\\w{1,3}$")) # extracts filename extension
  device(filename, width = width, height = height)
  print(plot)
  dev.off()
}
