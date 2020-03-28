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
library(raster)



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

ii <- sample(nrow(pr), 1000)

covs_clean <- covs_clean[ii, ]
pr <- pr[ii, ]
#+ trans

covs_clean <-
  covs_clean %>% 
  mutate(accessibility = log1p(accessibility),
         CHIRPS = log1p(CHIRPS),
         VIIRS = log1p(VIIRS))

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
                     max.edge = c(1, 20), 
                     cutoff = 1, 
                     min.angle = 21, 
                     offset = c(0.1, 30))

mesh$n

spde <- (inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	

coords <- 
  pr %>% filter(random_holdout == 0) %>% select(longitude, latitude) %>% as.matrix
Apix <- INLA::inla.mesh.project(mesh, loc = coords)$A



#+ setup_tmb, ache = FALSE

compile('spatially_varying_logit.cpp')

dyn.load(dynlib("spatially_varying_logit"))


parameters <- list(intercept = -5,
                   slope = rep(0, ncol(covs_clean)),
                   log_sigma = 0,
                   log_rho = 4,
                   nodemean = rep(0, nrow(spde$M0)),
                   log_covsigma = 0,
                   log_covrho = 4,
                   nodecov = rep(0, nrow(spde$M0))
                   
)

input_data <- list(x = as.matrix(covs_clean[pr$random_holdout == 0, ]),
                   Apixel = Apix,
                   spde = spde,
                   positive_cases = pr$pf_pos[pr$random_holdout == 0],
                   examined_cases = pr$examined[pr$random_holdout == 0],
                   priormean_intercept = -2,
                   priorsd_intercept = 2,
                   priormean_slope = 0,
                   priorsd_slope = 0.1,
                   prior_rho_min = 10,
                   prior_rho_prob = 0.01,
                   prior_sigma_max = 0.4,
                   prior_sigma_prob = 0.01,
                   prior_covrho_min = 20,
                   prior_covrho_prob = 0.01,
                   prior_covsigma_max = 0.1,
                   prior_covsigma_prob = 0.01,
                   nu = 1
                   
)

#+ fit_tmb

obj <- MakeADFun(
  data = input_data, 
  parameters = parameters,
  random = c('nodemean', 'nodecov'),
  DLL = "spatially_varying_logit")

its <- 1000
opt <- nlminb(obj$par, obj$fn, obj$gr, 
         control = list(iter.max = its, eval.max = 2*its, trace = 0))


#+ examined

pars <- split(obj$env$last.par.best, names(obj$env$last.par.best))

r <- raster('../../../data/covariates/housing/2019_Nature_Africa_Housing_2000.geotiff')
raster_pts <- rasterToPoints(r, spatial = TRUE)
coords_pred <- raster_pts@coords

Amatrix <- inla.mesh.project(mesh, loc = as.matrix(coords_pred))$A

field <- (Amatrix %*% pars$nodemean)[, 1]
field_ras <- rasterFromXYZ(cbind(coords_pred, field))
plot(field_ras)

covfield <- (Amatrix %*% pars$nodecov)[, 1]
covfield_ras <- rasterFromXYZ(cbind(coords_pred, covfield))
plot(covfield_ras)


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




