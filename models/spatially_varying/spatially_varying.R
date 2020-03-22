#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "Spatially varying coefficients model"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'spatially_varying'

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
library(TMB)



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
  mutate(month_start = pr$month_start) %>%
  mutate(year_start = pr$year_start) %>% 
  mutate(year_startrep = pr$year_start) %>% 
  mutate(country = as.numeric(factor(pr$country)))
  

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


outline.hull <- inla.nonconvex.hull(as.matrix(distinct(pr[, c('longitude', 'latitude')])), 
                                    convex = -0.02, 
                                    concave = -0.02,
                                    resolution = 400)
plot(outline.hull$loc, type = 'l')

#+ build_mesh

mesh <- inla.mesh.2d(pr[, c('longitude', 'latitude')], 
                     boundary = outline.hull,
                     max.edge = c(0.8, 20), 
                     cutoff = 0.8, 
                     min.angle = 21, 
                     offset = c(0.1, 30))

mesh$n




#+ setup_tmb

compile('spatially_varying_logit.cpp')

dyn.load(dynlib("spatially_varying_logit"))



#+ fit_tmb


#+ predict_base_random


summary_base_r <- summarise(pr$pf_pos[pr_inla$random_holdout == 1], 
                            pr_inla$examined[pr_inla$random_holdout == 1],
                            m1$summary.fitted$mean[1:nrow(pr_inla)][is.na(pr_inla$pf_pos)],
                            pr_inla[pr_inla$random_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0('base', name), 
                      covariates = 'base',
                      method = name,
                      cv = 'random',
                      mae = summary_base_r$weighted_mae,
                      correlation = summary_base_r$correlation,
                      time = m1$cpu.used[[2]])
summary_base_r$weighted_mae
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



ggplot(pr, aes(month_start, y = pf_pr, size = examined)) + 
  geom_point() +
  geom_smooth()

#' next bits

#+ ses

sessionInfo()




