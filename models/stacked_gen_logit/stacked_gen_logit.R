#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "stacked generalisation models with predictions transformed"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'stacked_gen_logit'

#' A fairly standard set of covariates
#' Monthly data.

#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, cache.lazy = FALSE,
                      fig.width = 8, fig.height = 5)

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

library(INLA)
library(INLAutils)
library(readr)



source('../../helper_functions/get_preds.R')
source('../../helper_functions/summarise.R')


if(!dir.exists('models')) dir.create('models')


#+ read_pr_data, eval = TRUE

pr <- fread("../../../data/derived/malariaAtlas_pr.csv")

#+ subset

pr <- pr %>% filter(year_start >= 2000)
pr <- pr %>% filter(continent == 'Africa')



#'# Base Data
#' ## Random CV


#+ read_in_models

dirs <- c('../nnet', 
          '../rf', 
          '../enet',
          '../ppr',
          '../xgboost')

models <- list()
weights <- rep(NA, length(dirs))
covs_clean <- list()

for(i in seq_along(dirs)){
  models[[i]] <- get(load(paste0(dirs[i], '/models/base_r.RData')))
  weights[i] <- min(models[[i]]$results$MAE, na.rm = TRUE)
  covs_clean[[i]] <- get_preds(models[[i]])
}

weights

weights <- (1 / weights) / sum(1 / weights)

weights_df <-
  data.frame(model = dirs %>% gsub('\\.\\.\\/', '', .), weight = weights)

kable(weights_df %>% arrange(desc(weights)), digits = 3)

covs_clean <- do.call(cbind, covs_clean) %>% as.data.frame




# data for final predictions
newdata <- list()
for(i in seq_along(dirs)){
  preds <- read_csv(paste0(dirs[i], '/random_base_errors.csv'))
  newdata[[i]] <- preds$pred
}
newdata <- do.call(cbind, newdata) %>% data.frame
names(newdata) <- names(covs_clean)




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


#+ combine
# combine covs_clean and new data

covs_clean2 <- matrix(NA, nrow = nrow(pr), ncol = ncol(covs_clean)) %>% data.frame

covs_clean2[pr$random_holdout == 1, ] <- newdata
covs_clean2[pr$random_holdout == 0, ] <- covs_clean
covs_clean2[covs_clean2 <= 0] <- 0.0001
covs_clean2 <- sapply(covs_clean2, qlogis) %>% data.frame
names(covs_clean2) <- names(covs_clean)


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
                                     cbind(b0 = 1, covs_clean2)))



#+ fit_mbg_random

h.spec <- list(theta=list(prior='pccor1', param=c(0, 0.9)))

fixed_names <- names(covs_clean %>% dplyr::select(-contains('start')))
fixed <- paste0('f(', fixed_names, ', model="clinear",range=c(0,Inf),initial=0)', collapse = ' + ')
form <- as.formula(paste(
  'pf_pos ~ b0 + 0 + f(field, model=spde, group = field.group, control.group = list(model="ar1", hyper=h.spec)) + ', 
  fixed))



m1 <- inla(form, data = inla.stack.data(stk.env), 
           family = 'binomial', 
           Ntrials = pr_inla$examined, 
           control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.env)),
           control.inla = list(int.strategy = 'eb', strategy = 'gaussian'), 
           num.threads = detectCores())


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





#' next bits

#+ ses

sessionInfo()




