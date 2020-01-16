#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'title: "bart models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'bart'

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

source('../../helper_functions/summarise.R')
source('../../helper_functions/extract_year.R')


#+ read_pr_data, eval = TRUE

pr <- fread("../../../data/derived/malariaAtlas_pr.csv")


#+ extract_covs, eval = FALSE


dirs <- c(
  '/home/tim/timz/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Monthly',  
  '/home/tim/timz/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_DiurnalDifference/5km/Monthly',
  '/home/tim/timz/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Night/5km/Monthly',
  '/home/tim/timz/mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/EVI_v6/5km/Monthly',
  '/home/tim/timz/mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCB_v6/5km/Monthly',
  '/home/tim/timz/mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCW_v6/5km/Monthly',
  '/home/tim/timz/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Monthly',
  '/home/tim/timz/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/5km/Monthly'
)
names1 <- c('LST_Day', 'LST_DiurnalDifference', 'LST_Night', 'EVI_v6',
            'TCB_v6', 'TCW_v6', 'TSI_Pf_Dynamic', 'CHIRPS'
)

static <- c(
  '/home/tim/timz/mastergrids/Other_Global_Covariates/Accessibility/Weiss/accessibility_to_cities_2015_v1.0.tif',
  '/home/tim/timz/mastergrids/Other_Global_Covariates/Aridity_v2/5km/Synoptic/Aridity_Index_v2.Synoptic.Overall.Data.5km.mean.tif',
  '/home/tim/timz/mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif',
  '/home/tim/timz/mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Annual_Clean_Composites/5km/Annual/VIIRS-DNB_Clean-Background.2016.Annual.Data.5km.mean.tif'
)
names2 <- c('accessibility', 'aridity', 'elevation', 'VIIRS')

covs <- rep(list(rep(NA, nrow(pr))), length(dirs))

for(i in seq_along(dirs)){
  
  
  short_files <- list.files(dirs[[i]], pattern = '\\.tif$')
  if(any(grepl('mean\\.5km\\.mean', short_files))){
    short_files <- grep('mean\\.5km\\.mean', short_files, value = TRUE)
  }
  year_vec <- str_extract(short_files, '20..\\.') %>% gsub('\\.', '', .) %>% as.numeric
  
  month_vec <- str_extract(short_files, '[0-9]{4}\\.[0-9]{2}') %>% gsub('[0-9]{4}\\.', '', .) %>% as.numeric
  
  
  full_files <-  list.files(dirs[[i]], pattern = '\\.tif$', full.names = TRUE)
  if(any(grepl('mean\\.5km\\.mean', full_files))){
    full_files <- grep('mean\\.5km\\.mean', full_files, value = TRUE)
  }
  
  if(month_vec[1] != 1){
    month_vec_old <- month_vec
    month_vec <- c(1:(month_vec[1] - 1), month_vec)
    year_vec <- c(rep(year_vec[1], (month_vec_old[1] - 1)), year_vec)
    full_files <- c(rep(full_files[1], (month_vec_old[1] - 1)), full_files)
  }
  
  covs[[i]] <- extract_year_month(coords = pr %>% dplyr::select(longitude, latitude),
                                  data_year = pr$year_start,
                                  data_month = pr$month_start,
                                  files = full_files,
                                  year_vec, 
                                  month_vec)  
}



covs_df <- do.call(cbind, covs) %>% data.frame
names(covs_df) <- names1

write.csv(covs_df, 'base_tmp.csv')



covs2 <- rep(list(rep(NA, nrow(pr))), length(static))
coords <- pr %>% dplyr::select(longitude, latitude) %>% as.matrix
for(i in seq_along(static)){
  r <- raster(static[i])
  covs2[[i]] <- raster::extract(r, coords)
}

covs_df2 <- do.call(cbind, covs2) %>% data.frame
names(covs_df2) <- names2




#+ clean_covs, eval = FALSE
covs <- cbind(covs_df, covs_df2)
sapply(covs, function(x) mean(is.na(x)))

covs_clean <- impute(covs)

#+ add_table_covs, eval = FALSE

covs_clean <- cbind(covs_clean, dplyr::select(pr, year_start, month_start))

#+ combine, eval = FALSE
covs_clean <- scale(covs_clean)
covs_clean <- as.data.frame(covs_clean)

#+ write_covs, eval = FALSE

write.csv(covs_clean, '../../../data/extracted_covs/base.csv')


if(!dir.exists('models')) dir.create('models')

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


#'# Base Data
#' ## Random CV


#+ fit_base_random, cache = TRUE

tg <- expand.grid(mtry = c(3, 8, 14), 
                  min.node.size = c(5, 20, 50),  
                  splitrule = 'variance')


m_base_r <- train(y = pr$pf_pr[pr$random_holdout == 0],
                x = covs_clean[pr$random_holdout == 0, ],
                method = 'ranger', 
                weights = pr$examined[pr$random_holdout == 0],
                tuneGrid = tg,
                metric = 'MAE',
                trControl = trainControl(method = 'boot632', number = 1)
)

save(m_base_r, file = 'models/base_r.RData')


#+ summary_base_random
plot(m_base_r)



#+ predict_base_random
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

write.csv(summary_base_r$errors, 'random_base_errors.csv')



#' next bits

#+ ses

sessionInfo()




