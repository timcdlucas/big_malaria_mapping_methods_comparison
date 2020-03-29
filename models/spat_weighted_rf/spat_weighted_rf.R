#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "spatially weighted rf models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'spat_weighted_rf'

#' A fairly standard set of covariates
#' Monthly data.

#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, fig.width = 8, fig.height = 5, cache.lazy = FALSE)

set.seed(110120)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tidyr)

library(ranger)

library(imputeMissings)
library(stringr)
library(knitr)

library(fields)



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




#'# Base Data
#' ## Random CV



#+ inner_cv, warning = FALSE

inner_holdout <- sample(which(pr$random_holdout == 0), 20)

partable <- expand.grid(N = c(1000, 4000), mtry = c(5, 8, 12))
mae <- rep(NA, length(nrow(partable)))

for(n in seq_len(nrow(partable))){
  
  print(n)
  N <- partable$N[n]
  preds <- rep(NA, inner_holdout %>% length)
  
  for(i in seq_along(inner_holdout)){
    
    #print(i)
    distances <- 
      fields.rdist.near(
        pr %>% slice(inner_holdout) %>% select(latitude, longitude)  %>% slice(i) %>% as.matrix,
        pr %>% filter(random_holdout == 0) %>% select(latitude, longitude)  %>% as.matrix,
        delta = 100,
        max.points = 5e5)
    
    this_indices <- order(distances$ra)[seq(N)]
    
    this_pr <- pr[pr$random_holdout == 0, ][this_indices, ]
    this_covs_clean <- covs_clean[pr$random_holdout == 0, ][this_indices, ]
    
    this_pr$weights <- this_pr$examined
    
    
    d <- data.frame(pf_pr = this_pr$pf_pr, this_covs_clean)
    this_model <- ranger(pf_pr ~ ., data = d , 
                         case.weights = this_pr$weights, 
                         mtry = partable$mtry[n],
                         min.node.size = 5)
    
    preds[i] <- predict(this_model, data = covs_clean[pr$random_holdout == 1, ][i, ], type = 'response')$predictions
  }
  
  mae[n] <- mean(abs(preds - pr$pf_pr[inner_holdout]))
}


ggplot(cbind(partable, mae = mae), aes(x = mtry, y = mae, colour = N)) + 
  geom_point()







#+ fit_base_random, cache = TRUE, results = 'hide', message = FALSE, warning = FALSE

N <- partable[which.min(mae), ]$N

geo_weight <- FALSE
preds <- rep(NA, pr$pf_pos[pr$random_holdout == 1] %>% length)



  for(i in seq_along(pr$pf_pos[pr$random_holdout == 1])){
    
    if((i %% 10000) == 0) print(i)

    distances <- 
      fields.rdist.near(
        pr %>% filter(random_holdout == 1) %>% select(latitude, longitude)  %>% slice(i) %>% as.matrix,
        pr %>% filter(random_holdout == 0) %>% select(latitude, longitude)  %>% as.matrix,
        delta = 100,
        max.points = 5e5)
    
    this_indices <- order(distances$ra)[seq(N)]
    
    this_pr <- pr[pr$random_holdout == 0, ][this_indices, ]
    this_covs_clean <- covs_clean[pr$random_holdout == 0, ][this_indices, ]
    
    if(geo_weight){
      this_pr$weights <- this_pr$examined * exp(- distances$ra[this_indices])
    } else {
      this_pr$weights <- this_pr$examined
    }
    
    d <- data.frame(pf_pr = this_pr$pf_pr, this_covs_clean)
    this_model <- ranger(pf_pr ~ ., data = d , 
                         case.weights = this_pr$weights,
                         min.node.size = 5,
                         mtry = partable[which.min(mae), ]$mtry)
    
    preds[i] <- predict(this_model, data = covs_clean[pr$random_holdout == 1, ][i, ], type = 'response')$predictions
  }




#+ predict_base_random, cache = FALSE
pred_base_r <- preds


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
                      time = NA)

write.csv(summary, 'random_base_summary.csv')
summary_base_r$weighted_mae

errors <- data.frame(name = paste0('base', name), 
                     covariates = 'base',
                     method = name,
                     pred = pred_base_r,
                     errors = summary_base_r$errors)

write.csv(errors, 'random_base_errors.csv')




#' next bits

#+ ses

sessionInfo()




