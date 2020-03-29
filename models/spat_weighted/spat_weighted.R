#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "nnet models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'spat_weighted'

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

inner_holdout <- sample(which(pr$random_holdout == 0), 500)
Nvec <- c(200, 500, 1000, 2000, 4000)
mae <- rep(NA, length(Nvec))

for(n in seq_along(Nvec)){
  
  N <- Nvec[n]
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
    this_model <- glm(pf_pr ~ ., data = d , weights = this_pr$weights, family = binomial())
    
    preds[i] <- predict(this_model, newdata = covs_clean[pr$random_holdout == 1, ][i, ], type = 'response')
  }
  
  mae[n] <- mean(abs(preds - pr$pf_pr[inner_holdout]))
}


plot(mae ~ Nvec)







#+ fit_base_random, cache = TRUE, results = 'hide', message = FALSE, warning = FALSE

N <- Nvec[which.min(mae)]
geo_weight <- FALSE
preds <- rep(NA, pr$pf_pos[pr$random_holdout == 1] %>% length)
coefficients <- data.frame(matrix(NA, ncol = ncol(covs_clean) + 1,
                                  nrow = pr$pf_pos[pr$random_holdout == 1] %>% length))

system.time(
for(i in seq_along(pr$pf_pos[pr$random_holdout == 1])){

  #print(i)
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
  this_model <- glm(pf_pr ~ ., data = d , weights = this_pr$weights, family = binomial())
  
  preds[i] <- predict(this_model, newdata = covs_clean[pr$random_holdout == 1, ][i, ], type = 'response')
  coefficients[i, ] <- coefficients(this_model)
}
)


#+ plots

coefficients_extra <- cbind(pr[pr$random_holdout == 1, ] %>% select(latitude, longitude), 
                            coefficients)
names(coefficients_extra) <- c('latitude', 'longitude', names(coefficients(this_model)))

coefficients_extra <- coefficients_extra %>% 
                        filter(!is.na(VIIRS))

coef_long <- pivot_longer(coefficients_extra,
                          cols = 3:ncol(coefficients_extra),
                          names_to = 'covariate',
                          values_to = 'estimate')

ggplot(coef_long, aes(longitude, latitude, colour = estimate)) +
  geom_point() + 
  facet_wrap(~ covariate) +
  scale_colour_viridis_c()



coefficients_extra <- cbind(pr[pr$random_holdout == 1, ] %>% select(latitude, longitude), 
                            scale(coefficients))
names(coefficients_extra) <- c('latitude', 'longitude', names(coefficients(this_model)))

coefficients_extra <- coefficients_extra %>% 
  filter(!is.na(VIIRS))

coef_long <- pivot_longer(coefficients_extra,
                          cols = 3:ncol(coefficients_extra),
                          names_to = 'covariate',
                          values_to = 'estimate')

ggplot(coef_long, aes(longitude, latitude, colour = estimate)) +
  geom_point() + 
  facet_wrap(~ covariate) +
  scale_colour_viridis_c()



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




