#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'title: "Base set of covariates for extra response data"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---


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
library(lubridate)
library(tidyr)
library(readr)

source('helper_functions/summarise.R')
source('helper_functions/extract_year.R')


#+ read_pr_data, eval = TRUE

mosq <- read.csv('../data/other_response/mosquito_occurrence.csv')
eir <- read.csv('../data/other_response/Monthly_EIR_data_for_Africa.csv')
files <- list.files('../data/other_response/espen', full.names = TRUE)

dfs <- lapply(files, read_csv)

lf <- do.call(rbind, dfs)
lf$Latitude <- as.numeric(lf$Latitude)
lf$Longitude <- as.numeric(lf$Longitude)

# Sort EIR

eir2 <-
  eir %>%
    pivot_longer(cols = starts_with('value'), names_to = 'month_column', values_to = 'eir' )

eir2 <-
  eir2 %>% 
    mutate(month_column_num = month_column %>% gsub('value', '', .) %>% as.numeric,
           month_start = Start.Month + month_column_num - 1) %>% 
    filter(month_start <= 12)

d <- data.frame(response = c(rep('mosq', nrow(mosq)), rep('eir', nrow(eir2)), rep('lf', nrow(lf))),
                latitude = c(mosq$latitude, eir2$`Latitude..degrees.`, lf$Latitude),
                longitude = c(mosq$longitude, eir2$`Longitude..degrees.`, lf$Longitude),
                year_start = as.numeric(c(mosq$year_start, eir2$Start.Year, lf$Year_start)),
                month_start = c(mosq$month_start, eir2$month_start, rep(6, nrow(lf))))

write.csv(d, '../data/other_response/combined_other_response.csv')

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

covs <- rep(list(rep(NA, nrow(d))), length(dirs))

for(i in seq_along(dirs)){
  print(i)
  
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
  
  covs[[i]] <- extract_year_month(coords = d %>% dplyr::select(longitude, latitude),
                                  data_year = d$year_start,
                                  data_month = d$month_start,
                                  files = full_files,
                                  year_vec, 
                                  month_vec)  
}



covs_df <- do.call(cbind, covs) %>% data.frame
names(covs_df) <- names1

write.csv(covs_df, 'base_extra_tmp.csv')



covs2 <- rep(list(rep(NA, nrow(d))), length(static))
coords <- d %>% dplyr::select(longitude, latitude) %>% as.matrix
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

covs_clean <- cbind(covs_clean, dplyr::select(d, year_start, month_start))

#+ combine, eval = FALSE
covs_clean <- scale(covs_clean)
covs_clean <- as.data.frame(covs_clean)

#+ write_covs, eval = FALSE


write.csv(covs_clean, 'extracted_covs/base_extra_response.csv')


#write.csv(covs_clean, '../../../data/extracted_covs/base_extra_response.csv')

