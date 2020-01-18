#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "mbg models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'mbg'

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
library(knitr)

library (INLA)
library(INLAutils)



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


if(Sys.info()["sysname"] != 'Windows'){
  message('using INLA unix workaround')
  INLA:::inla.dynload.workaround()
} else {
  message('Not using INLA unix workaround. Expect you are using winows.')
}


#+ setup_mbg_random

load('../../../data/derived/mesh.RData')


pr_inla <- pr
pr_inla$pf_pos[pr_inla$random_holdout == 1] <- NA
pr_inla$year_group <- as.numeric(cut(pr_inla$year_start, c(1999, 2005, 2010, 2015, 2030)))




# Define penalised complexity priors for random field. 
spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, 
                            prior.range = c(10, 0.01), prior.sigma = c(0.4, 0.01))
field.indices <- inla.spde.make.index("field", n.mesh = mesh$n, n.group = length(unique(pr_inla$year_group)))
Aest <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pr_inla[, c('longitude', 'latitude')]),
                         group = pr_inla$year_group)

stk.env <- inla.stack(tag = 'estimation', ## tag
                      data = list(pf_pos = pr_inla$pf_pos, examined = pr_inla$examined),
                      A = list(Aest, 1),  ## Projector matrix for space, fixed.
                      effects = list(field = field.indices,
                                     cbind(b0 = 1, covs_clean)))




#+ fit_mbg_random

fixed <- paste(names(covs_clean %>% dplyr::select(-contains('start'))), collapse = ' + ')
h.spec <- list(theta=list(prior='pccor1', param=c(0, 0.9)))
form <- as.formula(paste('pf_pos ~ b0 + 0 + f(field, model=spde, group = field.group, control.group = list(model="ar1", hyper=h.spec)) + ', fixed))



m1 <- inla(form, data = inla.stack.data(stk.env), 
           family = 'binomial', 
           Ntrials = pr_inla$examined, 
           control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.env)),
           control.inla = list(int.strategy = 'eb', strategy = 'gaussian'))


save(m1, file = 'models/inla1.RData')


#+ predict_base_random


summary_mbg_r <- summarise(pr$pf_pos[pr_inla$random_holdout == 1], 
                           pr_inla$examined[pr_inla$random_holdout == 1],
                           m1$summary.fitted$mean[1:nrow(pr_inla)][is.na(pr_inla$pf_pos)],
                           pr_inla[pr_inla$random_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0('base', name), 
                      covariates = 'base',
                      method = name,
                      cv = 'random',
                      mae = summary_mbg_r$weighted_mae,
                      correlation = summary_mbg_r$correlation,
                      time = m1$cpu.used[[2]])

write.csv(summary, 'random_mbg_summary.csv')


errors <- data.frame(name = paste0('base', name), 
                     covariates = 'base',
                     method = name,
                     pred = m1$summary.fitted$mean[1:nrow(pr_inla)][is.na(pr_inla$pf_pos)],
                     errors = summary_base_r$errors)

write.csv(errors, 'random_base_errors.csv')








#+ summary_base_random, cache = FALSE


autoplot(m1)

kable(m1$summary.fixed)

kable(m1$summary.hyper)

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







#' next bits

#+ ses

sessionInfo()




