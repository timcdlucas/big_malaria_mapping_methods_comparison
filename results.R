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

knitr::opts_chunk$set(cache = FALSE, fig.width = 8, fig.height = 5)

set.seed(110120)

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(knitr)
library(forcats)

#'# Covariates

#+ read_summaries, warning = FALSE, results = 'hide', message = FALSE

files <- list.files('covariates', recursive = TRUE, pattern = '_summary.csv', full.names = TRUE)

d <- lapply(files, read_csv)

data <- do.call(rbind, d)


files_big <- list.files('covariates', recursive = TRUE, pattern = '_errors.csv', full.names = TRUE)

d_big <- lapply(files, read_csv)

data_big <- do.call(rbind, d)



#+ plots

data %>% 
  filter(cv == 'random') %>% 
  ggplot(aes(x = fct_reorder(covariates, mae, .fun = min), 
                 fill = method, y = mae)) + 
    geom_bar(stat = 'identity', position = 'dodge')


data %>% 
  filter(cv == 'random') %>% 
  ggplot(aes(x = method, fill = covariates, y = mae)) + 
  geom_bar(stat = 'identity', position = 'dodge')



#+ table

data %>% 
  filter(cv == 'random') %>% 
  dplyr::select(-X1, -name, -cv, -time) %>% 
  kable(caption = 'Table of results for varying covariates',
        digits = 3)



data %>% 
  filter(cv == 'random') %>% 
  arrange(method) %>% 
  dplyr::select(-X1, -name, -cv, -time) %>% 
  kable(caption = 'Table of results for varying covariates',
        digits = 3)


#+ spatial, eval = FALSE
data %>% 
  filter(cv == 'spatial') %>% 
  ggplot(aes(x = covariates, fill = method, y = mae)) + 
    geom_bar(stat = 'identity', position = 'dodge')



#'# Models



#+ read_model_summaries, warning = FALSE, results = 'hide', message = FALSE

files <- list.files('models', recursive = TRUE, pattern = '_summary.csv', full.names = TRUE)

d <- lapply(files, read_csv)

data_mod <- do.call(rbind, d)

data_mod <- data %>% 
              filter(covariates == 'base') %>% 
              rbind(data_mod)


#+ read_error


files_big <- list.files('models', recursive = TRUE, pattern = '_errors.csv', full.names = TRUE)

d_big <- lapply(files_big, read_csv)

data_big <- do.call(rbind, d)


#+ plots_models

data_mod %>% 
  filter(cv == 'random') %>% 
  ggplot(aes(x = fct_reorder(method, mae), y = mae)) + 
    geom_bar(stat = 'identity', position = 'dodge')




#+ table_mod




data_mod %>% 
  filter(cv == 'random') %>% 
  arrange(mae) %>% 
  dplyr::select(-X1, -name, -cv, -time) %>% 
  kable(caption = 'Table of results for varying methods',
        digits = 3)



