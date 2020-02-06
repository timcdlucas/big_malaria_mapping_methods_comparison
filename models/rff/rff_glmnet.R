#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "GPs with random fourier features"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'rff_glmnet'

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
library(parallel)
library(doParallel)
library(knitr)
library(tictoc)
library(caret)

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
         VIIRS = log1p(VIIRS)) %>% 
  mutate(latitude = scale(pr$latitude), longitude = scale(pr$longitude))



#'# Base Data
#' ## Random CV


#+ fit_base_random, cache = TRUE, results = 'hide', message = FALSE

m <- list()
bw_vec <- c(0.3, 0.6, 0.8, 1, 1.2)

Y <- as.matrix(pr$pf_pos[pr$random_holdout == 0])
examined <- as.matrix(pr$examined[pr$random_holdout == 0])
X <- as.matrix(covs_clean)
K <- 1000
ncovs = ncol(X)


omega <- t(matrix(rnorm(K * ncovs), ncol = ncovs))
cl <- makeForkCluster(6)
registerDoParallel(cl)


for(i in seq_along(bw_vec)){
  
  feat <- X %*% (omega * bw_vec[i])
  
  rff <- sqrt(1 / K) * cbind(sin(feat), cos(feat))
  colnames(rff) <- paste('X', seq_len(ncol(rff)))
  
  
  
  m[[i]] <- train(y = pr$pf_pr[pr$random_holdout == 0],
                    x = rff[pr$random_holdout == 0, ],
                    method = 'glmnet', 
                    weights = pr$examined[pr$random_holdout == 0],
                    tuneLength = 15,
                    metric = 'MAE',
                    trControl = trainControl(method = 'boot', number = 1, 
                                             search = 'grid',
                                             savePredictions = TRUE)
  )
  
  plot(m_base_r)
  
}
  
maes <- sapply(m, function(x) min(x$results$MAE))

m_base_r <- m[[which.min(maes)]]
  
save(m_base_r, file = 'models/base_r_glmnet.RData')


stopCluster(cl)

#+ summary_base_random, cache = FALSE

plot(m_base_r)

plot(maes ~ bw_vec)
# Band width

#+ predict_base_random, cache = FALSE

feat <- X %*% (omega * bw_vec[which.min(maes)])

rff <- sqrt(1 / K) * cbind(sin(feat), cos(feat))
colnames(rff) <- paste('X', seq_len(ncol(rff)))

pred_base_r <- predict(m_base_r, newdata = rff[pr$random_holdout == 1, ])

pred_base_r[pred_base_r < 0] <- 0

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




