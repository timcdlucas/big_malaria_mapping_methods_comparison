#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'title: "Main data setup"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---


#+ setup
knitr::opts_chunk$set(cache = FALSE, fig.width = 8, fig.height = 5)
set.seed(1122)

#+ get_pr_data, eval = TRUE

pr <- getPR(country = 'all', species = 'pf')

pr <- 
  pr %>% 
    filter(!is.na(examined),
           !is.na(positive),
           !is.na(latitude),
           !is.na(longitude))


#+ create_random_holdout

random_holdout <- c(rep(0, floor(nrow(pr) * 0.8)), rep(1, ceiling(nrow(pr) * 0.2))) %>% sample



#+ create_spatial_holdout

spatial_folds <- kmeans(pr[, c('longitude', 'latitude')], 12)$cluster

WorldData <- map_data('world') %>% filter(region != "Antarctica") %>% fortify

cbind(pr, cluster = factor(spatial_folds)) %>% 
  ggplot(aes(x = longitude, y = latitude, colour = cluster)) + 
  geom_map(data = WorldData, map = WorldData,
             aes(x = long, y = lat, group = group, map_id=region),
             fill = NA, colour = "#7f7f7f", size=0.5) + 
    geom_point() +
    xlim(-80, 180) + 
    ylim(-35, 40) + 
    coord_equal() + 
    theme_minimal()


spatial_holdout <- spatial_folds %in% c(6, 4, 1)
cbind(pr, cluster = factor(spatial_holdout)) %>% 
  ggplot(aes(x = longitude, y = latitude, colour = cluster)) + 
  geom_map(data = WorldData, map = WorldData,
             aes(x = long, y = lat, group = group, map_id=region),
             fill = NA, colour = "#7f7f7f", size=0.5) + 
    geom_point() +
    xlim(-80, 180) + 
    ylim(-35, 40) + 
    coord_equal() + 
    theme_minimal()

spatial_holdout %>% sum
(!spatial_holdout) %>% sum

#+ Combine_and_write

pr <- cbind(pr, 
            random_holdout = random_holdout, 
            spatial_holdout = spatial_holdout)

write.csv(pr, "../data/derived/malariaAtlas_pr.csv")





