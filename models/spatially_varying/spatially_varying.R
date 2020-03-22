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

pr_inla <- pr
pr_inla$pf_pos[pr_inla$random_holdout == 1] <- NA
pr_inla$year_group <- as.numeric(cut(pr_inla$year_start, c(1999, 2005, 2010, 2015, 2030)))

# Define penalised complexity priors for random field. 
spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, 
                            prior.range = c(10, 0.01), prior.sigma = c(0.4, 0.01))

spde_cov <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, 
                            prior.range = c(10, 0.01), prior.sigma = c(0.1, 0.01))


field.indices <- inla.spde.make.index("field", n.mesh = mesh$n, n.group = length(unique(pr_inla$year_group)))

#covar1.indices <- inla.spde.make.index("cov1", n.mesh = mesh$n, n.group = length(unique(pr_inla$year_group)))

Ncovs <- 2
cov_indices_list <-
  lapply(seq_len(Ncovs), function(i) 
    inla.spde.make.index(paste0("cov", i), n.spde = mesh$n))
    

Aest <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pr_inla[, c('longitude', 'latitude')]),
                         group = pr_inla$year_group)

Acov_list <- 
  lapply(seq_len(Ncovs), function(i)
    inla.spde.make.A(mesh = mesh, loc = as.matrix(pr_inla[, c('longitude', 'latitude')]),
                     weights = covs_clean[, i]))
#Ac1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pr_inla[, c('longitude', 'latitude')]),
                         #group = pr_inla$year_group, weights = covs_clean$LST_Day)

Alist <- c(Aest, Acov_list, rep(1, ncol(covs_clean) + 1))
effectslist <- c(field = list(field.indices),
                 cov1 = cov_indices_list,
                 cbind(b0 = 1, covs_clean))

#effectslist <- list(field = field.indices,
#                    cov1 = covar1.indices,
#                    cbind(b0 = 1, covs_clean))
stk.env <- inla.stack(tag = 'estimation', ## tag
                      data = list(pf_pos = pr_inla$pf_pos, examined = pr_inla$examined),
                      A = Alist,  ## Projector matrix for space, fixed.
                      effects = effectslist)



#+ fit_mbg_random

h.spec <- list(theta=list(prior='pccor1', param=c(0, 0.9)))
cov_ar1 <- list(theta=list(prior='pccor1', param=c(0, 0.95)))

hyper.rw2 <- list(prec = list(prior="pc.prec", param = c(1, 0.01)))
hyper.country <- list(prec = list(prior="pc.prec", param = c(0.4, 0.01)))


fixed <- paste(names(covs_clean %>% dplyr::select(-contains('start'))), collapse = ' + ')

form1 <- 'pf_pos ~ b0 + 0 + '
form2 <- 'f(field, model = spde, group = field.group, control.group = list(model="ar1", hyper=h.spec)) + '
form3 <- 'f(year_start, model="rw2", hyper = hyper.rw2, scale.model = TRUE) + '
#form4 <- 'f(year_startrep, model="rw2", hyper = hyper.country, scale.model = TRUE, replicate = country) + '
formc1 <- 'f(cov1, model = spde_cov, group = cov1.group, control.group = list(model="ar1", hyper=cov_ar1)) + '
form_spat_var <- lapply(seq_len(Ncovs), function(i)
  paste0('f(cov', i, 
         ', model = spde_cov, group = cov', i, 
         '.group, control.group = list(model="ar1", hyper=cov_ar1)) + '))
form_spat_var <- paste0(form_spat_var, collapse = '')

form <- as.formula(paste(form1, form2, form3, form_spat_var, fixed))
form <- as.formula(paste(form1, form3, form_spat_var, fixed))



m1 <- inla(form, data = inla.stack.data(stk.env), 
           family = 'binomial', 
           Ntrials = pr_inla$examined, 
           control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.env)),
           control.inla = list(int.strategy = 'eb', strategy = 'gaussian', adjust.weights = FALSE),
           num.threads = 8, verbose = T)


save(m1, file = 'models/inla1.RData')


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




