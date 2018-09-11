knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(message = T)
knitr::opts_chunk$set(warning = F)
knitr::opts_chunk$set(error = T)
knitr::opts_chunk$set(cache = F)
Sys.setlocale("LC_ALL", "en_IE.UTF-8")

ipak(tidyverse)
ipak(mlr)
ipak(parallelMap)
ipak(parallel)
# newly added 2018-09-05:
ipak(plotmo)
ipak(reshape2)
ipak(gplots)
ipak(gridExtra)
ipak(plotly)
ipak(ComplexHeatmap)
ipak(circlize)

#Custom functions
walk(list.files("./R", pattern = "*.R", full.names = T), source)

# dynamically assign number of processors
# if (system('hostname', intern = TRUE) == 'spinoza.embl.de') {
#    parallelStartMulticore(cpus = 10)
# } else {
#    parallelStartMulticore(cpus = detectCores() - 1)
# }
