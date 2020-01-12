#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'title: "Results"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---


#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, fig.width = 8, fig.height = 5)

set.seed(110120)

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

#'# Covariates

#+ read_summaries

files <- list.files('covariates', recursive = TRUE, pattern = '_summary.csv', full.names = TRUE)

d <- lapply(files, read_csv)

data <- do.call(rbind, d)


#+ plots

data %>% 
  filter(cv == 'random') %>% 
  ggplot(aes(x = covariates, fill = method, y = mae)) + 
    geom_bar(stat = 'identity', position = 'dodge')




#'# Models




