#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "bart models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'bart'

#' A fairly standard set of covariates
#' Monthly data.

#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, fig.width = 8, fig.height = 5)

set.seed(110120)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(imputeMissings)
library(stringr)
library(parallel)
library(caret)
library(doParallel)

source('../../helper_functions/summarise.R')
source('../../helper_functions/extract_year.R')

if(!dir.exists('models')) dir.create('models')


#+ read_pr_data, eval = TRUE

pr <- fread("../../../data/derived/malariaAtlas_pr.csv")



#+ read_back_in

covs_clean <- fread('../../../data/extracted_covs/base.csv')
covs_clean <- covs_clean %>% dplyr::select(-V1)


#+ subset
covs_clean <- covs_clean[pr$year_start >= 2000, ]
pr <- pr %>% filter(year_start >= 2000)


covs_clean <- covs_clean[pr$continent == 'Africa', ]
pr <- pr %>% filter(continent == 'Africa')


#+ trans

covs_clean <-
  covs_clean %>% 
  mutate(accessibility = log1p(accessibility),
         CHIRPS = log1p(CHIRPS),
         VIIRS = log1p(VIIRS))


#+ setup_parallel, cache= FALSE




#'# Base Data bart
#' ## Random CV


#+ fit_base_random, cache = TRUE, results = 'hide', message = FALSE

options( java.parameters = "-Xmx48g" )

N <- 20
gr <- data.frame(num_trees = seq(30, 60, length.out = N), 
                  k = runif(N, min = 0, max = 5),
                  alpha = runif(N, min = .9, max = 1),
                  beta = runif(N, min = 0, max = 4),
                  nu = runif(N, min = 0, max = 5))


m_base_r <- train(y = pr$pf_pr[pr$random_holdout == 0],
                x = covs_clean[pr$random_holdout == 0, ],
                method = 'bartMachine', 
                weights = pr$examined[pr$random_holdout == 0],
                tuneGrid = gr,
                metric = 'MAE',
                trControl = trainControl(method = 'cv', number = 3, 
                                         search = 'random', 
                                         selectionFunction = 'oneSE'),
                mem_cache_for_speed = FALSE,
                run_in_sample = FALSE,
                serialize = TRUE
)

save(m_base_r, file = 'models/base_r.RData')





#+ summary_base_random, cache = FALSE

plot(m_base_r$results$MAE ~ m_base_r$results$num_trees)
plot(m_base_r$results$MAE ~ m_base_r$results$beta)
plot(m_base_r$results$MAE ~ m_base_r$results$nu)
plot(m_base_r$results$MAE ~ m_base_r$results$alpha)
plot(m_base_r$results$MAE ~ m_base_r$results$k)

kable(m_base_r$results[, 1:8], digits = 2)

#plot(m_base_r)



#+ predict_base_random, cache = FALSE
pred_base_r <- predict(m_base_r, newdata = covs_clean[pr$random_holdout == 1, ])


summary_base_r <- summarise(pr$pf_pos[pr$random_holdout == 1], 
                          pr$examined[pr$random_holdout == 1],
                          pred_base_r,
                          pr[pr$random_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0('base', name), 
                      covariates = 'base',
                      method = name,
                      cv = 'random',
                      mae = summary_base_r$weighted_mae,
                      correlation = summary_base_r$correlation,
                      time = m_base_r$times$everything[[1]])

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




