#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'title: "Add world bank covariates"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'world_bank'

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
#library(binomTools)

source('../../helper_functions/summarise.R')
source('../../helper_functions/extract_year.R')


#+ read_pr_data, eval = TRUE

pr <- fread("../../../data/derived/malariaAtlas_pr.csv")

#+ base

covs_clean <- fread('../../../data/extracted_covs/base.csv')
covs_clean <- covs_clean %>% dplyr::select(-V1)



#+ extract_covs

files <- list.files('../../../data/covariates/WorldBankIndicators_Input/', 
                    pattern = '\\.csv$',
                    full.names = TRUE)

# Don't want rural and urban pop percent.
files <- files[1:(length(files) - 1)]

raw <- lapply(files, function(x) fread(x, skip = 2, header = TRUE))

rawlong <- lapply(raw, function(x) pivot_longer(x, cols = -c(1:4), names_to = 'year'))

data <- do.call(rbind, rawlong)
data <- data %>% mutate(year = as.numeric(year))

data$transform <- 'positive'
percent <- 
  data$`Indicator Name` %in% 
    c("Immunization, DPT (% of children ages 12-23 months)",
      "Access to electricity (% of population)",
      "Health expenditure, public (% of total health expenditure)",
      "Health expenditure, total (% of GDP)",
      "Literacy rate, adult total (% of people ages 15 and above)",
      "Pregnant women receiving prenatal care (%)",
      "Primary completion rate, total (% of relevant age group)",
      "Rural population (% of total population)"
      )

identity <-   
  data$`Indicator Name` %in% 
    c("GDP growth (annual %)",
      "GDP per capita (current US$)")
   
data$transform[percent] <- 'percent'
data$transform[identity] <- 'identity'
data$transform %>% table

# Make transformed vars

  
data <- 
  data %>% 
    group_by(`Country Name`, `Indicator Code`) %>% 
    mutate(all_non_na = all(is.na(value))) %>%
    filter(!all_non_na) %>% 
    mutate(interpolate = spline(value, x = year, xout = year)$y) %>% 
    ungroup

data <-
  data %>% 
    mutate(interpolate =
             case_when(
               transform == 'positive' & interpolate < 0 ~ 0,
               transform == 'percent' & interpolate < 0 ~ 0,
               transform == 'percent' & interpolate > 100 ~ 100,
               TRUE ~ interpolate
             ))



data_wide <-
  data %>% 
    rename(year_column = year) %>%
    pivot_wider(names_from = `Indicator Code`, values_from = interpolate, - c(`Indicator Name`, value, transform)) 


#+ join_back

wrld_covs <- 
  pr %>% 
    left_join(data_wide, by = c('country_id' = 'Country Code', 'year_start' = 'year_column')) %>% 
    dplyr::select(SH.MED.CMHW.P3:SP.RUR.TOTL.ZS)


#+ clean_covs
sapply(wrld_covs, function(x) mean(is.na(x)))

wrld_covs <- wrld_covs[, -1]

wrld_covs <- impute(wrld_covs)

#+ combine, eval = FALSE
wrld_covs <- scale(wrld_covs)
wrld_covs <- as.data.frame(wrld_covs)

#+ write_covs, eval = FALSE

write.csv(covs_clean, '../../../data/extracted_covs/world_bank.csv')


#+ bind_covs, eval = FALSE
covs_clean <- cbind(covs_clean, wrld_covs)



#+ fit_models

if(!dir.exists('models')) dir.create('models')


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


#+ basic_plots

covs_clean %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram() +
  facet_wrap(~ name, scale = 'free')




#' # Elastic net
#' ## Random CV


#+ fit_enet_random, cache = TRUE

cl <- makeForkCluster(8)
registerDoParallel(cl)

m_enet_r <- train(y = pr$pf_pr[pr$random_holdout == 0],
                  x = covs_clean[pr$random_holdout == 0, ],
                  method = 'enet', 
                  weights = pr$examined[pr$random_holdout == 0],
                  tuneLength = 10,
                  metric = 'MAE',
                  trControl = trainControl(method = 'boot632', number = 1)
)
stopCluster(cl)

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
                     pred = pred_base_r,
                     errors = summary_base_r$errors)

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
                     pred = pred_base_r,
                     errors = summary_base_r$errors)

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
                     pred = pred_base_r,
                     errors = summary_base_r$errors)

write.csv(errors, 'random_mbg_errors.csv')



#' next bits

#+ ses

sessionInfo()




