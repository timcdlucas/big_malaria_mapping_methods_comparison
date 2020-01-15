#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'title: "Build global mesh"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, fig.width = 8, fig.height = 5)

#+ libs

library(INLA)
library(data.table)
library(dplyr)
library(raster)
library(INLAutils)
library(ggplot2)

#+ data

pr <- fread("../data/derived/malariaAtlas_pr.csv")
pr <- pr %>% filter(continent == 'Africa')

#+ build_outline



outline.hull <- inla.nonconvex.hull(as.matrix(distinct(pr[, c('longitude', 'latitude')])), 
                                    convex = -0.02, 
                                    concave = -0.02,
                                    resolution = 400)
plot(outline.hull$loc, type = 'l')

#+ build_mesh

mesh <- inla.mesh.2d(pr[, c('longitude', 'latitude')], 
                     boundary = outline.hull,
                     max.edge = c(0.6, 20), 
                     cutoff = 0.6, 
                     min.angle = 21, 
                     offset = c(0.1, 30))

mesh$n






#+ plot_mesh

autoplot(mesh)


WorldData <- map_data('world') %>% filter(region != "Antarctica") %>% fortify

autoplot(mesh) +
  geom_map(data = WorldData, map = WorldData,
             aes(x = long, y = lat, group = group, map_id=region),
             fill = NA, colour = "red", size=0.5) + 
    xlim(-30, 60) + 
    ylim(-50, 50) + 
    coord_equal() + 
    theme_minimal()


#+ save
save(mesh, file = '../data/derived/mesh.RData')


