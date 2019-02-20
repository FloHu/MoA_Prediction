knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(message = T)
knitr::opts_chunk$set(warning = F)
knitr::opts_chunk$set(error = T)
knitr::opts_chunk$set(cache = F)
Sys.setlocale("LC_ALL", "en_IE.UTF-8")

ipak(magrittr)
ipak(tidyverse)
ipak(mlr)
ipak(parallelMap)
ipak(parallel)
ipak(plotmo)
ipak(reshape2)
ipak(gplots)
ipak(gridExtra)
ipak(plotly)
ipak(ComplexHeatmap)
ipak(circlize)
ipak(googledrive)
ipak(viridis)

#Custom functions
walk(list.files("./R", pattern = "*.R", full.names = T), source)

# change plot defaults:
theme_set(theme_bw())
theme_update(text = element_text(size = 12))

