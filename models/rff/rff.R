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

name <- 'rff'

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
library(greta)
library(doParallel)
library(knitr)
library(tictoc)

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
  mutate(latitude = pr$latitude, longitude = pr$longitude)



#'# Base Data
#' ## Random CV


#+ fit_base_random, cache = TRUE, results = 'hide', message = FALSE

tic()
Y <- as.matrix(pr$pf_pos[pr$random_holdout == 0])
examined <- as.matrix(pr$examined[pr$random_holdout == 0])
X <- as.matrix(covs_clean[pr$random_holdout == 0, ])
K <- 2500
ncovs = ncol(X)


# priors
tau <- exponential(1)
bw <- gamma(1,1)
sigma <- gamma(1,1)
W <- normal(0, 1, dim = c(2 * K, 1)) * tau


# Do the math
omega <- t(matrix(rnorm(K * ncovs), ncol = ncovs))

feat <- X %*% (omega * bw)

rff <- sqrt(1 / K) * cbind(sin(feat), cos(feat))


mu <- rff %*% W
prev <- ilogit(mu)
distribution(Y) <- binomial(examined, prev)

model <- model(W, tau, bw)

m_base_r <- opt(model, optimiser = adam())
time <- toc()

save(m_base_r, file = 'models/base_r.RData')


#+ summary_base_random, cache = FALSE

# Band width
m_base_r$par$bw

# Scale I think
1 / m_base_r$par$bw

Xpred <- expand.grid(rep(list(seq(-3, 3, 0.1)), 2)) %>% 
           cbind(matrix(0, nrow = nrow(.), ncol = ncovs - 2)) %>% 
            as.matrix
featpred <- Xpred %*% (omega*bw)
rffpred <- sqrt(1/K)*cbind(sin(featpred), cos(featpred))
mupred <- rffpred %*% W
prevpred <- ilogit(mupred)

Ypred <- calculate(prevpred, m_base_r$par)

data.frame(Xpred, Y = Ypred) %>%
  ggplot(aes(Var1, Var2, fill = Y)) + 
  geom_tile()



Xpred <- expand.grid(rep(list(seq(-3, 3, 0.1)), 2)) %>% 
  cbind(matrix(0, nrow = nrow(.), ncol = ncovs - 2), .) %>% 
  as.matrix
featpred <- Xpred %*% (omega*bw)
rffpred <- sqrt(1/K)*cbind(sin(featpred), cos(featpred))

mupred <- rffpred %*% W
prevpred <- ilogit(mupred)

Ypred <- calculate(prevpred, m_base_r$par)

data.frame(Xpred, Y = Ypred) %>%
  ggplot(aes(Var1, Var2, fill = Y)) + 
  geom_tile()

Ypred_insample <- calculate(prev, m_base_r$par)
preds <- data.frame(obs = pr$pf_pr[pr$random_holdout == 0],
                    pred = Ypred_insample)
ggplot(preds, aes(obs, pred)) + 
  geom_point() + 
  geom_smooth() + 
  geom_abline(slope = 1, intercept = 0)










#+ predict_base_random, cache = FALSE
Xpred <- covs_clean[pr$random_holdout == 1, ]
featpred <- Xpred %*% (omega*bw)
rffpred <- sqrt(1/K)*cbind(sin(featpred), cos(featpred))
mupred <- rffpred %*% W
prevpred <- ilogit(mupred)
pred_base_r <- calculate(prevpred, m_base_r$par)[, 1, drop = TRUE]


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
                      time = time$toc - time$tic)

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




