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
library(broom)
library(sjstats)

#+ read_in_pr
pr <- fread("../data/derived/malariaAtlas_pr.csv")
pr <- pr %>% filter(year_start >= 2000)
pr <- pr %>% filter(continent == 'Africa')

#'# Covariates

#+ read_summaries, warning = FALSE, results = 'hide', message = FALSE

files <- list.files('covariates', recursive = TRUE, pattern = '_summary.csv', full.names = TRUE)

d <- lapply(files, read_csv)

data <- do.call(rbind, d)


#files_big <- list.files('covariates', recursive = TRUE, pattern = '_errors.csv', full.names = TRUE)#

#d_big <- lapply(files, read_csv)

#data_big <- do.call(rbind, d)



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

data_mod <- data_mod %>% 
              filter(covariates == 'base')


#+ read_error


files_big <- list.files('models', recursive = TRUE, pattern = '_errors.csv', full.names = TRUE)

d_big <- lapply(files_big, read_csv)

data_big <- do.call(rbind, d_big)
data_big <- cbind(data_big, examined = pr$examined[pr$random_holdout == 1])



#+ plots_models

data_mod %>% 
  filter(cv == 'random') %>% 
  ggplot(aes(x = fct_reorder(method, mae), y = mae)) + 
    geom_bar(stat = 'identity', position = 'dodge')




#+ table_mod

best_model <- data_mod$method[which.min(data_mod$mae)]
best_errors <- data_big$errors[data_big$method == best_model]
data_mod$p_value <- NA
method_vec <- unique(data_mod$method)

for(i in seq_along(method_vec)){
  if(method_vec[i] == best_model){
    data_mod$p_value[i] <- NA
  } else {
    d <- data_big %>% filter(method %in% c(method_vec[i], best_model))
    data_mod$p_value[i] <- wtd_mwu(d, errors, method, examined)$p.value
  }
}



data_mod %>% 
  filter(cv == 'random') %>% 
  arrange(mae) %>% 
  dplyr::select(method, mae, p_value, correlation) %>% 
  kable(caption = 'Table of results for varying methods',
        digits = 3)



