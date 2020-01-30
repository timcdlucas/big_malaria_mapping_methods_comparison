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
library(rstan)
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

m3 <- stan_model("rff3.stan", auto_write = TRUE)
m4 <- stan_model("rff4.stan", auto_write = TRUE)


Nmodels <- 36
par_table <- data.frame(k = sample(50:200, Nmodels), bw = exp(runif(Nmodels, -4, 0.2)), l1 = runif(Nmodels))

mae_inner <- rep(NA, Nmodels)
par(mfrow = c(sqrt(Nmodels), sqrt(Nmodels)))
for(i in seq_along(mae_inner)){
  
  lm1 <- lm.fit(as.matrix(covs_clean[, 1:12]), pr$pf_pr)
  data = list(y = lm1$residuals,
              x = covs_clean,
              n = nrow(pr),
              m = ncol(covs_clean),
              l1 = par_table$l1[i],
              bw = 0.1,
              k = par_table$k[i],
              omega = matrix(rnorm(par_table$k[i] * ncol(covs_clean)), ncol = ncol(covs_clean)))
  
  fitmap = optimizing(m4, data = data)
  ypredgp <- fitmap$par[grepl('fhat', names(fitmap$par))]
  
  ypred <- ypredgp + lm1$fitted.values
  
  mae_inner[i] <- weighted.mean(abs(ypred - pr$pf_pr), weight = pr$examined)
  title <- paste0('k ', par_table$k[i], ', bw ', round(par_table$bw[i], 2), ', l1 ',  round(par_table$l1[i], 2), ', mae ', round(mae_inner[i], 2))
  plot(pr$pf_pr, ypred, main = title)
  
}

ggplot(cbind(par_table, mae = mae_inner), aes(x = k, y = mae, colour = bw, size = l1)) + 
  geom_point() + 
  geom_smooth()

ggplot(cbind(par_table, mae = mae_inner), aes(x = bw, y = mae, colour = l1, size = k)) + 
  geom_point() + 
  geom_smooth()


ggplot(cbind(par_table, mae = mae_inner), aes(x = l1, y = mae, colour = k, size = bw)) + 
  geom_point() + 
  geom_smooth()


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




