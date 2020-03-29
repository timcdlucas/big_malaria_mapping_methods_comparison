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

name <- 'lf_rf'

#' A fairly standard set of covariates
#' Monthly data.

#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, fig.width = 8, fig.height = 5, cache.lazy = FALSE)

set.seed(110120)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(imputeMissings)
library(stringr)
library(parallel)
library(caret)
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

covs_clean$year_start <- covs_clean$year_start - median(covs_clean$year_start)

pr <- 
  pr %>% 
  select(latitude, longitude, examined, pf_pos, pf_pr, random_holdout)

lf <- lf %>% 
  mutate(random_holdout = 0) %>% 
  rename(latitude = Latitude, longitude = Longitude, 
         examined = Examined, pf_pos = Positive, pf_pr = Prevalence) %>% 
  select(latitude, longitude, examined, pf_pos, pf_pr, random_holdout)
pr <- bind_rows(pr, lf)

#'# Base Data bart
#' ## Random CV


#+ fit_base_random, cache = TRUE, results = 'hide', message = FALSE


N <- 20


m <- list()

wts <- rep(c(0.6, 0.8, 1), 2)
N <- 12
tuneGrid <- data.frame(mtry = sample(14, N, replace = TRUE), 
                       splitrule = sample(c("variance", "extratrees", "maxstat"), N, replace = TRUE), 
                       min.node.size = sample(40, N, replace = TRUE))

split <- 'malaria'

nmalaria <- nrow(covs_clean[pr$random_holdout == 0 & covs_clean$malaria == 1, ])
holdout <- sample(nmalaria, 7000)
train <- which(!(seq_len(nrow(covs_clean[pr$random_holdout == 0, ])) %in% holdout))

trCtrl<- 
  trainControl(indexOut = list(holdout),
               index = list(train),
               allowParallel = FALSE,
               savePredictions = TRUE)

for(i in seq_along(wts)){
  
  split <- if(i > (length(wts) / 2)) split <- NULL
  weight_vec <- pr$examined
  weight_vec[covs_clean$malaria == 0] <- weight_vec[covs_clean$malaria == 0] * wts[i]
  
  m[[i]] <- train(y = pr$pf_pr[pr$random_holdout == 0],
                  x = covs_clean[pr$random_holdout == 0, ],
                  method = 'ranger', 
                  weights = weight_vec[pr$random_holdout == 0],
                  tuneGrid = tuneGrid,
                  metric = 'MAE',
                  trControl = trCtrl,
                  always.split.variables = 'malaria')
}

perf <- sapply(m, function(x) x$results$MAE %>% min)


m_base_r <- m[[which.min(perf)]]
save(m_base_r, file = 'models/base_r.RData')





#+ summary_base_random, cache = FALSE

plot(perf)

plot(m_base_r$results$MAE ~ m_base_r$results$mtry)
plot(m_base_r$results$MAE ~ m_base_r$results$min.node.size)

kable(m_base_r$results[, seq_len(length(m_base_r$bestTune) + 3)], digits = 2)

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


summary_base_r$weighted_mae
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




