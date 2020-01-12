


summarise <- function(positive, examined, pred, coords){

  # Weighted Mae

  pr <- positive / examined
  error <- pr - pred
  abs_error <- abs(error)

  weighted_mae <- weighted.mean(abs_error, examined)


  # correlation

  correlation <- cor(pr, pred)

  # scatter

  p <- 
    data.frame(obs = pr, pred = pred, weight = examined) %>% 
      ggplot(aes(y = pred, x = pr, size = examined)) + 
        geom_point(alpha = 0.5) + 
        geom_smooth(show.legend = FALSE) + 
        geom_abline(slope = 1, intercept = 0)

  print(p)


  # spatial residuals


  # map of predictions 
  # not a raster because makes function very complicated.

  WorldData <- map_data('world') %>% filter(region != "Antarctica") %>% fortify

  p2 <- 
    cbind(coords, error) %>% 
      ggplot(aes(x = longitude, y = latitude, colour = error)) + 
      geom_map(data = WorldData, map = WorldData,
                 aes(x = long, y = lat, group = group, map_id=region),
                 fill = NA, colour = "#7f7f7f", size=0.5) + 
        geom_point() +
        xlim(-80, 180) + 
        scale_colour_distiller(palette = 'Spectral') +
        ylim(-35, 40) + 
        coord_equal() + 
        theme_minimal()
   

  # Combine outputs

  outs <- list(weighted_mae = weighted_mae, 
               correlation = correlation,
               errors = abs_error)

  return(outs)

}
