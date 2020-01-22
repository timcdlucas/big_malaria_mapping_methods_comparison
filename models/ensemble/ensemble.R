#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "Ensemble models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'ensemble'

#' A fairly standard set of covariates
#' Monthly data.

#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, cache.lazy = FALSE,
                      fig.width = 8, fig.height = 5)

set.seed(110120)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(readr)

library(imputeMissings)
library(stringr)
library(parallel)
library(caret)
library(doParallel)
library(knitr)
library(penalized)


source('../../helper_functions/summarise.R')
source('../../helper_functions/extract_year.R')
source('../../helper_functions/get_preds.R')

if(!dir.exists('models')) dir.create('models')


#+ read_pr_data, eval = TRUE

pr <- fread("../../../data/derived/malariaAtlas_pr.csv")

#+ subset

pr <- pr %>% filter(year_start >= 2000)
pr <- pr %>% filter(continent == 'Africa')



#'# Base Data
#' ## Random CV


#+ read_in_models

dirs <- c('../nnet', 
          '../rf', 
          '../enet',
          '../ppr',
          '../xgboost')

models <- list()
weights <- rep(NA, length(dirs))
covs_clean <- list()

for(i in seq_along(dirs)){
  models[[i]] <- get(load(paste0(dirs[i], '/models/base_r.RData')))
  weights[i] <- min(models[[i]]$results$MAE, na.rm = TRUE)
  covs_clean[[i]] <- get_preds(models[[i]])
}

weights

weights <- (1 / weights) / sum(1 / weights)

weights_df <-
  data.frame(model = dirs %>% gsub('\\.\\.\\/', '', .), weight = weights)

kable(weights_df %>% arrange(desc(weights)), digits = 3)

covs_clean <- do.call(cbind, covs_clean) %>% as.data.frame




# data for final predictions
newdata <- list()
for(i in seq_along(dirs)){
  preds <- read_csv(paste0(dirs[i], '/random_base_errors.csv'))
  newdata[[i]] <- preds$pred
}
newdata <- do.call(cbind, newdata) %>% data.frame
names(newdata) <- names(covs_clean)



#+ predict_base_random, cache = FALSE
#m_base_r <- lm(pf_pr ~ 0 + V1 + V2 + V3 + V4 + V5,
#              data = cbind(covs_clean,
#                           pf_pr = pr$pf_pr[pr$random_holdout == 0]))

m_base_r <- penalized(pf_pr, ~ V1 + V2 + V3 + V4 + V5, ~ 0,
                      positive = TRUE,
                      data = cbind(covs_clean, 
                                   pf_pr = pr$pf_pr[pr$random_holdout == 0]))


data.frame(model = dirs %>% gsub('\\.\\.\\/', '', .),
           weight = m_base_r@penalized) %>% 
  kable(digits = 3)



pred_base_r <- predict(m_base_r, penalized = newdata)[, 1]
#pred_base_r <- predict(m_base_r, newdata = newdata)

summary_base_r <- summarise(pr$pf_pos[pr$random_holdout == 1], 
                          pr$examined[pr$random_holdout == 1],
                          pred_base_r,
                          pr[pr$random_holdout == 1, c('longitude', 'latitude')])
summary_base_r$weighted_mae


summary <- data.frame(name = paste0('base', name), 
                      covariates = 'base',
                      method = name,
                      cv = 'random',
                      mae = summary_base_r$weighted_mae,
                      unweighted_mae = summary_base_r$unweighted_mae,
                      correlation = summary_base_r$correlation,
                      time = NA)


write.csv(summary, 'random_base_summary.csv')


errors <- data.frame(name = paste0('base', name), 
                     covariates = 'base',
                     method = name,
                     pred = pred_base_r,
                     errors = summary_base_r$errors)

write.csv(errors, 'random_base_errors.csv')


#' next bits

#+ ses

sessionInfo()




