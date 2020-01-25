#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "Spatial rf models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'enet_spatial'

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
library(caret)
library(doParallel)
library(fields)
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
         VIIRS = log1p(VIIRS))


#+ make_distances

coords <- pr %>% dplyr::select(longitude, latitude)
kmeans <- kmeans(coords, 500)

WorldData <- map_data('world') %>% filter(region != "Antarctica") %>% fortify

kmeans$centers %>% 
  as.data.frame %>% 
  ggplot(aes(x = longitude, y = latitude)) + 
  geom_map(data = WorldData, map = WorldData,
           aes(x = long, y = lat, group = group, map_id=region),
           fill = NA, colour = "#7f7f7f", size=0.5) + 
  geom_point() +
  lims(x = range(pr$longitude), y = range(pr$latitude)) +
  coord_equal()


dist_mat <- apply(coords, 1, function(x) 
  sqrt((x[[1]] - kmeans$centers[, 1]) ^ 2 + (x[[2]] - kmeans$centers[, 2]) ^ 2))

min_val <- min(dist_mat[dist_mat > 0])
dist_mat[dist_mat == 0] <- min_val / 2


covs_clean <- cbind(covs_clean, 1 / t(dist_mat))


#'# Base Data 
#' ## Random CV


#+ fit_base_random, cache = TRUE, results = 'hide', message = FALSE



cl <- makeForkCluster(10)
registerDoParallel(cl)

N <- 10


m_base_r <- train(y = pr$pf_pr[pr$random_holdout == 0],
                x = covs_clean[pr$random_holdout == 0, ],
                method = 'enet', 
                weights = pr$examined[pr$random_holdout == 0],
                tuneLength = N,
                metric = 'MAE',
                trControl = trainControl(method = 'cv', number = 3, 
                                         search = 'grid', 
                                         selectionFunction = 'oneSE'))

save(m_base_r, file = 'models/base_r.RData')


stopCluster(cl)



#+ summary_base_random, cache = FALSE

plot(m_base_r)

kable(m_base_r$results[, seq_len(length(m_base_r$bestTune) + 3)], digits = 2)



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


errors <- data.frame(name = paste0('base', name), 
                     covariates = 'base',
                     method = name,
                     pred = pred_base_r,
                     errors = summary_base_r$errors)

write.csv(errors, 'random_base_errors.csv')




#' next bits

#+ ses

sessionInfo()




