#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'title: "First test set of covariates and models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'weiss'

#' The covariates from Weiss 2015
#' Monthly data.
#' Transforms and interactions already made and then covariate selection performed.

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

source('../../helper_functions/summarise.R')
source('../../helper_functions/extract_year.R')


#+ read_pr_data, eval = TRUE

pr <- fread("../../../data/derived/malariaAtlas_pr.csv")


#+ extract_covs, eval = FALSE

dirs <- list.dirs('~/timz/mastergrids/Other_Global_Covariates/Pf_Covariates/')
dirs <- grep('Monthly$', dirs, value = TRUE)

covs <- rep(list(rep(NA, nrow(pr))), length(dirs))
for(i in seq_along(dirs)){
  
  print(i)
  
  short_files <- list.files(dirs[[i]], pattern = '\\.tif$')
  year_vec <- str_extract(short_files, '20..\\.') %>% gsub('\\.', '', .) %>% as.numeric
  month_vec <- str_extract(short_files, '..\\.Data') %>% gsub('\\.Data', '', .) %>% as.numeric
  
  full_files <-  list.files(dirs[[i]], pattern = '\\.tif$', full.names = TRUE)
  
  covs[[i]] <- extract_year_month(coords = pr %>% dplyr::select(longitude, latitude),
                                  data_year = pr$year_start,
                                  data_month = pr$month_start,
                                  files = full_files,
                                  year_vec, 
                                  month_vec)
  
}


covs_df <- do.call(cbind, covs) %>% data.frame
names(covs_df) <- str_extract(dirs, 'PF_V[0-9]+\\/') %>% gsub('\\/', '', .)


#write.csv(covs_df, 'weiss_tmp.csv')



#+ clean_covs, eval = FALSE
covs <- covs_df
sapply(covs, function(x) mean(is.na(x)))

covs_clean <- impute(covs)

#+ add_table_covs, eval = FALSE

covs_clean <- cbind(covs_clean, dplyr::select(pr, year_start, month_start))

#+ combine, eval = FALSE
covs_clean <- scale(covs_clean)
covs_clean <- as.data.frame(covs_clean)

#+ write_covs, eval = FALSE

write.csv(covs_clean, '../../../data/extracted_covs/weiss.csv')


#+ read_back_in

covs_clean <- fread('../../../data/extracted_covs/weiss.csv')
covs_clean <- covs_clean %>% dplyr::select(-V1)

#+ basic_plots

covs_clean %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = value)) + 
    geom_histogram() +
    facet_wrap(~ name, scale = 'free')




#' # Elastic net
#' ## Random CV


#+ fit_enet_random

m_enet_r <- train(y = pr$pf_pr[pr$random_holdout == 0],
                  x = covs_clean[pr$random_holdout == 0, ],
                  method = 'enet', 
                  weights = pr$examined[pr$random_holdout == 0],
                  tuneLength = 10,
                  metric = 'MAE',
                  trControl = trainControl(method = 'boot632', number = 3, 
                                           selectionFunction = 'oneSE')
                  )


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

write.csv(summary_enet_r$errors, 'random_enet_errors.csv')


#' ## Spatial CV

#+ fit_enet_spatial

m_enet_s <- train(y = pr$pf_pr[pr$spatial_holdout == 0],
                  x = covs_clean[pr$spatial_holdout == 0, ],
                  method = 'enet', 
                  weights = pr$examined[pr$spatial_holdout == 0],
                  tuneLength = 10,
                  metric = 'MAE',
                  trControl = trainControl(method = 'boot632', number = 3, 
                                           selectionFunction = 'oneSE')
)

#+ predict_enet_spatial

pred_enet_s <- predict(m_enet_s, newdata = covs_clean[pr$spatial_holdout == 1, ])


summary_enet_s <- summarise(pr$pf_pos[pr$spatial_holdout == 1], 
                            pr$examined[pr$spatial_holdout == 1],
                            pred_enet_s,
                            pr[pr$spatial_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0(name, 'enet'), 
                      covariates = name,
                      method = 'enet',
                      cv = 'spatial',
                      mae = summary_enet_s$weighted_mae,
                      correlation = summary_enet_s$correlation,
                      time = m_enet_s$times$everything[[1]])

write.csv(summary, 'spatial_enet_summary.csv')

write.csv(summary_enet_s$errors, 'spatial_enet_errors.csv')


#'# Random forest
#' ## Random CV


#+ fit_rf_random, eval = TRUE

tg <- expand.grid(mtry = c(7, 12, 17), 
                  min.node.size = c(5, 20, 50),  
                  splitrule = 'variance')


m_rf_r <- train(y = pr$pf_pr[pr$random_holdout == 0],
                  x = covs_clean[pr$random_holdout == 0, ],
                  method = 'ranger', 
                  weights = pr$examined[pr$random_holdout == 0],
                  tuneGrid = tg,
                  metric = 'MAE',
                  trControl = trainControl(method = 'boot632', number = 3, 
                                           selectionFunction = 'oneSE')
                  )


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

write.csv(summary_rf_r$errors, 'random_rf_errors.csv')



#' ## Spatial CV


#+ fit_rf_spatial, eval = TRUE

m_rf_s <- train(y = pr$pf_pr[pr$spatial_holdout == 0],
                x = covs_clean[pr$spatial_holdout == 0, ],
                method = 'ranger', 
                weights = pr$examined[pr$spatial_holdout == 0],
                tuneGrid = tg,
                metric = 'MAE',
                trControl = trainControl(method = 'boot632', number = 3, 
                                         selectionFunction = 'oneSE', allowParallel = FALSE)
)


#+ summary_rf_spatial

plot(m_rf_s)


#+ predict_rf_spatial

pred_rf_s <- predict(m_rf_s, newdata = covs_clean[pr$spatial_holdout == 1, ])

summary_rf_s <- summarise(pr$pf_pos[pr$spatial_holdout == 1], 
                          pr$examined[pr$spatial_holdout == 1],
                          pred_rf_s,
                          pr[pr$spatial_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0(name, 'rf'), 
                      covariates = name,
                      method = 'rf',
                      cv = 'spatial',
                      mae = summary_rf_s$weighted_mae,
                      correlation = summary_rf_s$correlation,
                      time = m_rf_s$times$everything[[1]])

write.csv(summary, 'spatial_rf_summary.csv')

write.csv(summary_rf_s$errors, 'spatial_rf_errors.csv')




#'# mbg
#' ## Random CV


#+ setup_mbg_random

load('../../../data/derived/mesh.RData')

Aest <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pr[, c('longitude', 'latitude')]))

# Define penalised complexity priors for random field. 
spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, prior.range = c(10, 0.01), prior.sigma = c(0.4, 0.01))


stk.env <- inla.stack(tag = 'estimation', ## tag
                      data = list(pf_pos = pr$pf_pos, examined = pr$examined),
                      A = list(Aest, 1),  ## Projector matrix for space, fixed, then sum to one fixed.
                      effects = list(space = seq_len(spde$n.spde), 
                                     cbind(b0 = 1, covs_clean)))



#+ fit_mbg_random
	
fixed <- paste(names(covs_clean %>% dplyr::select(-contains('start'))), collapse = ' + ')
form <- as.formula(paste('pf_pos ~ b0 + 0 + f(space, model = spde) + ', fixed))

m1 <- inla(form, data = inla.stack.data(stk.env), family = 'binomial', 
           Ntrials = examined, 
           control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.env)),
           num.threads = 4)

save(m1, file = 'models/inla1.RData')

#+ predict_mbg_random

m1$summary.random$space %>% 
  cbind(mesh$loc) %>% 
  ggplot(aes(`1`, `2`, colour = plogis(mean + m1$summary.fixed['b0', 'mean']))) + 
    geom_point() + 
    scale_colour_viridis_c()

#' ## Spatial CV

#+ fit_mbg_spatial


#+ predict_mbg_spatial


#'# Stacked generalisation
#' ## Random CV

#+ fit_stackgen_random


#+ predict_stackgen_random


#+ fit_stackgen_spatial


#' ## Spatial CV
#+ predict_stackgen_spatial





