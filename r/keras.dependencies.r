# Dependency Loader, Dien Wu, initially written by Ben Fasoli for STILTv2

####
if (!'smurf_wd' %in% ls()) smurf_wd <- getwd()

### Source R scripts 
rsc <- dir(file.path(smurf_wd, 'r/src'), pattern = '.*\\.r$',
           full.names = T, recursive = T)

rsc <- c(rsc, dir(file.path(smurf_wd, 'r/src'), pattern = '.*\\.R$',
              full.names = T, recursive = T))
              
invisible(lapply(rsc, source))

### Load external libraries
if (!'lib.loc' %in% ls()) lib.loc <- NULL
libs <- load_libs('dplyr', 'parallel', 'raster', 'ggplot2', 'ggmap', 
                  'ggpubr', 'reshape2', 'stringr', 'rasterVis', 'lutz',
                  'rslurm', 'lubridate', 'gridExtra', lib.loc = lib.loc)

# load additional packages for Keras
library(reticulate); library(keras); library(tensorflow)    
