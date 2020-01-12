#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'title: "Feature  of covariates and models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'feature_expansion'

#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, fig.width = 8, fig.height = 5)

set.seed(110120)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(malariaAtlas)
library(raster)
library(imputeMissings)
library(caret)

source('../../helper_functions/summarise.R')


#+ read_pr_data, eval = TRUE

pr <- read.csv("../../../data/derived/malariaAtlas_pr.csv")


#+ read_covs

r_files <- list.files("../../../data/covariates/first_test/", pattern = '\\.tif$', full.name = TRUE)

r <- lapply(r_files, raster)

r <- r[-5]

sapply(r, extent)



#+ extract covs

vals <- rep(list(rep(NA, nrow(pr))), length(r))

for(i in seq_along(r)){
  print(i)
  vals[[i]] <- extract(r[[i]], pr[, c('longitude', 'latitude')])
}

covs <- do.call(cbind, vals) %>% data.frame
names(covs) <- sapply(r, names)


#+ clean_covs

sapply(covs, function(x) mean(is.na(x)))

covs_clean <- impute(covs)
covs_clean <- log1p(covs_clean)

#+ feature_engineer

covs_clean <- cbind(covs_clean, log(covs_clean - min(covs_clean) + 0.0001), covs_clean ^ 2)
names(covs_clean) <- c(names(covs), paste0('log_', names(covs)), paste0('square_', names(covs))) 

covs_clean <- scale(covs_clean)
covs_clean <- as.data.frame(covs_clean)




#+ basic_plots

covs_clean %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = value)) + 
    geom_histogram() +
    facet_wrap(~ name, scale = 'free')



#+ write_covs

write.csv(covs_clean, paste0('../../../data/extracted_covs/', name, '.csv'))

#+ create_random_holdout

random_holdout <- c(rep(0, floor(nrow(pr) * 0.8)), rep(1, ceiling(nrow(pr) * 0.2))) %>% sample


#+ create_spatial_holdout


#' Elastic net


#+ fit_enet_random

m_enet_r <- train(y = pr$pr[random_holdout == 0],
                  x = covs_clean[random_holdout == 0, ],
                  method = 'enet', 
                  weights = pr$examined[random_holdout == 0],
                  tuneLength = 10,
                  metric = 'MAE',
                  trControl = trainControl(method = 'cv', number = 5, 
                                           selectionFunction = 'oneSE')
                  )


#+ summary_enet_random

plot(m_enet_r)



#+ predict_enet_random

pred_enet_r <- predict(m_enet_r, newdata = covs_clean[random_holdout == 1, ])


summary_enet_r <- summarise(pr$positive[random_holdout == 1], 
                            pr$examined[random_holdout == 1],
                            pred_enet_r,
                            pr[random_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0(name, 'enet'), 
                      covariates = name,
                      method = 'enet',
                      cv = 'random',
                      mae = summary_enet_r$weighted_mae,
                      correlation = summary_enet_r$correlation,
                      time = m_enet_r$times$everything[[1]])

write.csv(summary, 'random_enet_summary.csv')

write.csv(summary_enet_r$errors, 'random_enet_errors.csv')



#+ fit_enet_spatial


#+ predict_enet_spatial



#'# Random forest


#+ fit_rf_random

tg <- expand.grid(mtry = c(3, 7, 14, 20), 
                  min.node.size = c(1, 5, 10, 20),  
                  splitrule = 'variance')


m_rf_r <- train(y = pr$pr[random_holdout == 0],
                  x = covs_clean[random_holdout == 0, ],
                  method = 'ranger', 
                  weights = pr$examined[random_holdout == 0],
                  tuneGrid = tg,
                  metric = 'MAE',
                  trControl = trainControl(method = 'cv', number = 5, 
                                           selectionFunction = 'oneSE')
                  )


#+ summary_rf_random

plot(m_rf_r)



#+ predict_rf_random

pred_rf_r <- predict(m_rf_r, newdata = covs_clean[random_holdout == 1, ])


summary_rf_r <- summarise(pr$positive[random_holdout == 1], 
                          pr$examined[random_holdout == 1],
                          pred_rf_r,
                          pr[random_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0(name, 'rf'), 
                      covariates = name,
                      method = 'rf',
                      cv = 'random',
                      mae = summary_rf_r$weighted_mae,
                      correlation = summary_rf_r$correlation,
                      time = m_rf_r$times$everything[[1]])

write.csv(summary, 'random_rf_summary.csv')

write.csv(summary_rf_r$errors, 'random_rf_errors.csv')




#+ fit_rf_spatial


#+ predict_rf_spatial



#'# mbg

#+ fit_mbg_random


#+ predict_mbg_random


#+ fit_mbg_spatial


#+ predict_mbg_spatial


#'# Stacked generalisation

#+ fit_stackgen_random


#+ predict_stackgen_random


#+ fit_stackgen_spatial


#+ predict_stackgen_spatial





