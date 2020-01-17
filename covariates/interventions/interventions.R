#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'title: "Add intervention covariates"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'intervention'

#' A fairly standard set of covariates
#' Monthly data.

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
library(stringr)
library(parallel)
library(caret)
library(INLA)
library(data.table)
library(raster)

source('../../helper_functions/summarise.R')
source('../../helper_functions/extract_year.R')


#+ read_pr_data, eval = TRUE

pr <- fread("../../../data/derived/malariaAtlas_pr.csv")

#+ base

covs_base <- fread('../../../data/extracted_covs/base.csv')
covs_base <- covs_base %>% dplyr::select(-V1)



#+ extract_covs

rs <- getRaster(c('Indoor residual spraying (IRS) coverage', 
                 'Insecticide-treated bednet (ITN) coverage'),
               year = list(2000:2015, 2000:2015))

files <- list.files('../../../data/covariates/housing/', pattern = '\\.geotiff$', full.names = TRUE)

house <- raster::stack(files)

coords <- pr %>% dplyr::select(longitude, latitude) %>% as.matrix

irs <- extract_year(coords,   
                    pr$year_start,
                    rs[[seq_along(2000:2015)]],
                    2000:2015)
itn <- extract_year(coords,   
                    pr$year_start,
                    rs[[seq_along(2000:2015) + 16]],
                    2000:2015)

housing <- extract_year(coords,   
                        pr$year_start,
                        house,
                        c(2000, 2015))

covs <- data.frame(itn = itn, irs = irs, housing = housing)




#+ clean_covs
sapply(covs, function(x) mean(is.na(x)))
covs_clean <- impute(covs)
sapply(covs_clean, function(x) mean(is.na(x)))

#+ write_covs, eval = TRUE

write.csv(covs, paste0('../../../data/extracted_covs/', name, '.tif'))


if(!dir.exists('models')) dir.create('models')

#+ add_table_covs

covs_clean <- cbind(covs_base, covs)
covs_clean <- impute(covs_clean)


#+ subset
covs_clean <- covs_clean[pr$year_start >= 2000, ]
pr <- pr %>% filter(year_start >= 2000)


covs_clean <- covs_clean[pr$continent == 'Africa', ]
pr <- pr %>% filter(continent == 'Africa')


#+ basic_plots

covs_clean %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram() +
  facet_wrap(~ name, scale = 'free')




#' # Elastic net
#' ## Random CV


#+ fit_enet_random, cache = TRUE

m_enet_r <- train(y = pr$pf_pr[pr$random_holdout == 0],
                  x = covs_clean[pr$random_holdout == 0, ],
                  method = 'enet', 
                  weights = pr$examined[pr$random_holdout == 0],
                  tuneLength = 10,
                  metric = 'MAE',
                  trControl = trainControl(method = 'boot632', number = 1)
)

save(m_enet_r, file = 'models/enet_r.RData')

#+ summary_enet_random

plot(m_enet_r)



#+ predict_enet_random

pred_enet_r <- predict(m_enet_r, newdata = covs_clean[pr$random_holdout == 1, ])


summary_enet_r <- summarise(pr$pf_pos[pr$random_holdout == 1], 
                            pr$examined[pr$random_holdout == 1],
                            pred_enet_r,
                            pr[pr$random_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0(name, 'enet'), 
                      covariates = name,
                      method = 'enet',
                      cv = 'random',
                      mae = summary_enet_r$weighted_mae,
                      correlation = summary_enet_r$correlation,
                      time = m_enet_r$times$everything[[1]])

write.csv(summary, 'random_enet_summary.csv')


errors <- data.frame(name = paste0(name, 'enet'), 
                     covariates = name,
                     method = 'enet',
                     pred = pred_enet_r,
                     errors = summary_enet_r$errors)

write.csv(errors, 'random_enet_errors.csv')



#'# Random forest
#' ## Random CV


#+ fit_rf_random, cache = TRUE

tg <- expand.grid(mtry = c(3, 8, 14), 
                  min.node.size = c(5, 20, 50),  
                  splitrule = 'variance')


m_rf_r <- train(y = pr$pf_pr[pr$random_holdout == 0],
                x = covs_clean[pr$random_holdout == 0, ],
                method = 'ranger', 
                weights = pr$examined[pr$random_holdout == 0],
                tuneGrid = tg,
                metric = 'MAE',
                trControl = trainControl(method = 'boot632', number = 1)
)

save(m_rf_r, file = 'models/rf_r.RData')


#+ summary_rf_random
plot(m_rf_r)



#+ predict_rf_random
pred_rf_r <- predict(m_rf_r, newdata = covs_clean[pr$random_holdout == 1, ])


summary_rf_r <- summarise(pr$pf_pos[pr$random_holdout == 1], 
                          pr$examined[pr$random_holdout == 1],
                          pred_rf_r,
                          pr[pr$random_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0(name, 'rf'), 
                      covariates = name,
                      method = 'rf',
                      cv = 'random',
                      mae = summary_rf_r$weighted_mae,
                      correlation = summary_rf_r$correlation,
                      time = m_rf_r$times$everything[[1]])

write.csv(summary, 'random_rf_summary.csv')



errors <- data.frame(name = paste0(name, 'rf'), 
                     covariates = name,
                     method = 'rf',
                     pred = pred_rf_r,
                     errors = summary_rf_r$errors)

write.csv(errors, 'random_rf_errors.csv')






#'# mbg
#' ## Random CV

#+ setup_inla

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


#+ predict_mbg_random


summary_mbg_r <- summarise(pr$pf_pos[pr_inla$random_holdout == 1], 
                           pr_inla$examined[pr_inla$random_holdout == 1],
                           m1$summary.fitted$mean[1:nrow(pr_inla)][is.na(pr_inla$pf_pos)],
                           pr_inla[pr_inla$random_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0(name, 'mbg'), 
                      covariates = name,
                      method = 'mbg',
                      cv = 'random',
                      mae = summary_mbg_r$weighted_mae,
                      correlation = summary_mbg_r$correlation,
                      time = m1$cpu.used[[2]])

write.csv(summary, 'random_mbg_summary.csv')



errors <- data.frame(name = paste0(name, 'mbg'), 
                     covariates = name,
                     method = 'mbg',
                     pred = m1$summary.fitted$mean[1:nrow(pr_inla)][is.na(pr_inla$pf_pos)],
                     errors = summary_mbg_r$errors)

write.csv(errors, 'random_mbg_errors.csv')



#' next bits

#+ ses

sessionInfo()




