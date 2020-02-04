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

name <- 'rff'

#' A fairly standard set of covariates
#' Monthly data.

#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, fig.width = 8, fig.height = 5)

set.seed(110120)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tidyr)

library(imputeMissings)
library(stringr)
library(parallel)
library(greta)
library(doParallel)
library(knitr)

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


#+ fit_base_random, cache = TRUE, results = 'hide', message = FALSE


Y <- as.matrix(pr$pf_pr[pr$random_holdout == 0])

X <- as.matrix(covs_clean[pr$random_holdout == 0, ])
K <- 500
ncovs = ncol(X)


omega <- t(matrix(rnorm(K * ncovs), ncol = ncovs))

tau <- exponential(1)
bw <- gamma(1,1)
sigma <- gamma(1,1)

W <- normal(0, 1, dim = c(2 * K, 1)) * tau
feat<-X %*% (omega * bw)

rff <- sqrt(1 / K) * cbind(sin(feat), cos(feat))
mu <- rff %*% W
distribution(Y) <- normal(mu, sigma)

model <- model(W, tau, bw)

MAP <- opt(model, optimiser = adam())


Xpred <- expand.grid(rep(list(seq(-3, 3, 0.1)), 2)) %>% 
           cbind(matrix(0, nrow = nrow(.), ncol = ncovs - 2)) %>% 
            as.matrix
featpred <- Xpred %*% (omega*bw)
rffpred <- sqrt(1/K)*cbind(sin(featpred), cos(featpred))
Ypred <- rffpred %*% W
Ypred <- calculate(Ypred, MAP$par)

data.frame(Xpred, Y = Ypred) %>%
  ggplot(aes(Var1, Var2, fill = Y)) + 
  geom_tile()



Xpred <- expand.grid(rep(list(seq(-3, 3, 0.1)), 2)) %>% 
  cbind(matrix(0, nrow = nrow(.), ncol = ncovs - 2), .) %>% 
  as.matrix
featpred <- Xpred %*% (omega*bw)
rffpred <- sqrt(1/K)*cbind(sin(featpred), cos(featpred))
Ypred <- rffpred %*% W
Ypred <- calculate(Ypred, MAP$par)

data.frame(Xpred, Y = Ypred) %>%
  ggplot(aes(Var1, Var2, fill = Y)) + 
  geom_tile()

Ypred_insample <- calculate(mu, MAP$par)
preds <- data.frame(obs = pr$pf_pr[pr$random_holdout == 0],
                    pred = Ypred_insample)
ggplot(preds, aes(obs, pred)) + 
  geom_point() + 
  geom_smooth() + 
  geom_abline(slope = 1, intercept = 0)



save(m_base_r, file = 'models/base_r.RData')




#+ summary_base_random, cache = FALSE


m_base_r$results %>% 
  dplyr::select(-RMSE, -Rsquared, -RMSESD, -MAESD, -RsquaredSD) %>% 
  pivot_longer(-MAE) %>% 
  ggplot(aes(x = value, y = MAE)) +
    geom_point() + 
    geom_smooth() + 
    facet_wrap(~ name, scale = 'free')

kable(arrange(m_base_r$results[, seq_len(length(m_base_r$bestTune) + 3)], desc(MAE)), digits = 2)

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




