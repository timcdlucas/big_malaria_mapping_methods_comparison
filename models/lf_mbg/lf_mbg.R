#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "rf models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'lf_mbg'

#' A fairly standard set of covariates
#' Monthly data.

#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, fig.width = 8, fig.height = 5)

set.seed(110120)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(imputeMissings)
library(stringr)
library(parallel)
library(INLA)
library(doParallel)
library(readr)
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
         VIIRS = log1p(VIIRS)) %>% 
  mutate(year_start = pr$year_start)



#+ get_lf

files <- list.files('../../../data/other_response/espen', full.names = TRUE)
stopifnot(length(files) > 0)

dfs <- lapply(files, read_csv)

lf <- do.call(rbind, dfs)
lf$Latitude <- as.numeric(lf$Latitude)
lf$Longitude <- as.numeric(lf$Longitude)
lf$Year_start[lf$Year_start == 'null'] <- NA
lf$Year_start <- as.character(lf$Year_start)

lf <- 
  lf %>% 
    mutate(Year_start = ifelse(Year_start %in% c(1998, 1999), 2000, Year_start)) %>% 
    mutate(Examined = as.numeric(Examined),
           Positive = as.numeric(Positive)) %>% 
    filter(!is.na(Year_start)) %>% 
    filter(Positive <= Examined) %>%
    filter(Year_start >= 2000) %>% 
    filter(!is.na(Examined), !is.na(Prevalence), !is.na(Positive)) %>% 
    filter(!is.na(Latitude), !is.na(Longitude))


dim(lf)



#+ lf_covs

library(fields)

distances <- 
  fields.rdist.near(
    lf %>% select(Latitude, Longitude) %>% as.matrix,
    pr %>% select(latitude, longitude)  %>% as.matrix,
    delta = 25,
    max.points = 5e8)

dd <- distances$ind %>% as.data.frame

dd$dist <- distances$ra

close <- 
  dd %>%
    group_by(V1) %>% 
    arrange(dist) %>% 
    slice(1)

rm(distances)
rm(dd)


lf_covs <- covs_clean[close$V2, ] %>% 
             dplyr::select(-month_start, -year_start) %>% 
             mutate(month_start = rep(0, nrow(close)),
                    year_start = as.numeric(lf$Year_start))

covs_clean$malaria <- 1
lf_covs$malaria <- 0

covs_clean <- bind_rows(covs_clean, lf_covs)



pr <- 
  pr %>% 
    select(latitude, longitude, examined, pf_pos, pf_pr, random_holdout, country, year_start)

lf <- lf %>% 
  mutate(random_holdout = 0,
         year_start = as.numeric(Year_start)) %>% 
  rename(latitude = Latitude, longitude = Longitude, 
         examined = Examined, pf_pos = Positive, pf_pr = Prevalence,
         country = Country) %>% 
  select(latitude, longitude, examined, pf_pos, pf_pr, random_holdout, year_start, country)
pr <- bind_rows(pr, lf)



covs_clean <-
  covs_clean %>% 
  mutate(month_start = pr$month_start) %>%
  mutate(year_start = pr$year_start) %>% 
  mutate(year_start_lf = ifelse(covs_clean$malaria == 0, pr$year_start, NA)) %>% 
  mutate(year_startrep = pr$year_start) %>% 
  mutate(country = as.numeric(factor(pr$country))) %>% 
  mutate(lf_iid = ifelse(malaria == 0, 1, NA))


#'# Base Data bart
#' ## Random CV



#+ fit_base_random, cache = FALSE, results = 'hide', message = FALSE


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
hyper.rw2 <- list(prec = list(prior="pc.prec", param = c(1, 0.01)))
hyper.country <- list(prec = list(prior="pc.prec", param = c(0.4, 0.01)))
hyper.iid <- list(prec = list(prior="pc.prec", param = c(1, 0.01)))

form1 <- 'pf_pos ~ b0 + 0 + '
form2 <- 'f(field, model=spde, group = field.group, control.group = list(model="ar1", hyper=h.spec)) + '
form3 <- 'f(year_start, model="rw2", hyper = hyper.rw2, scale.model = TRUE) + '
form4 <- 'f(year_startrep, model="rw2", hyper = hyper.country, scale.model = TRUE, replicate = country) + '
form5 <- 'f(year_start_lf, model="rw2", hyper = hyper.rw2, scale.model = TRUE) + '
form6 <- 'f(lf_iid, model = "iid", hyper = hyper.iid) + '

form <- as.formula(paste(form1, form2, form3, form4, form5, form6, fixed))



m1 <- inla(form, data = inla.stack.data(stk.env), 
           family = 'binomial', 
           Ntrials = pr_inla$examined, 
           control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.env)),
           control.inla = list(int.strategy = 'eb', strategy = 'gaussian'),
           num.threads = 8, verbose = TRUE)


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

summary_base_r$weighted_mae


errors <- data.frame(name = paste0('base', name), 
                     covariates = 'base',
                     method = name,
                     pred = m1$summary.fitted$mean[1:nrow(pr_inla)][is.na(pr_inla$pf_pos)],
                     errors = summary_base_r$errors)

write.csv(errors, 'random_base_errors.csv')







#+ summary_base_random, cache = FALSE


#autoplot(m1)

kable(m1$summary.fixed)

kable(m1$summary.hyper)


#' next bits

#+ ses

sessionInfo()




