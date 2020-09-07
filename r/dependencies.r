# Dependency Loader, Dien Wu, initially written by Ben Fasoli for STILTv2

####
if (!'smurf_wd' %in% ls()) smurf_wd <- getwd()

### Source R scripts 
rsc <- dir(file.path(smurf_wd, 'r/src'), pattern = '.*\\.r$',
           full.names = T, recursive = T)

rsc <- c(rsc, dir(file.path(smurf_wd, 'r/src'), pattern = '.*\\.R$',
              full.names = T, recursive = T))
rsc <- rsc[!grepl('keras', rsc)]    # remove strips depending on Keras
            
invisible(lapply(rsc, source))

### Load external libraries
if (!'lib.loc' %in% ls()) lib.loc <- NULL
libs <- load_libs('dplyr', 'parallel', 'raster', 'ggplot2', 'ggmap', 
                  'ggpubr', 'reshape2', 'stringr', 'rasterVis', 'lutz',
                  'rslurm', 'lubridate', 'gridExtra', 'rworldmap', 
                  lib.loc = lib.loc)

