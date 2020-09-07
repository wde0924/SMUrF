### grab TROPOMI SIF downscaled by Alex Turner, 
# written by DW, 10/24/2019 

# based on Alex's latest SIF file, need to convert SIF from vectors to grids
# return a data.frame

# 02/04/2020, DW, calculate the mean SIF of multiple TSIF files if nhrs != NULL 

grab.tsif <- function(tsif.path, timestr, minlon, maxlon, minlat, maxlat, 
                      nhrs = NULL, form = c('data.frame', 'raster')[2]) {

    library(ncdf4)

    # select the ones that need for scaling GPP 
    tsif.files <- list.files(path = tsif.path, pattern = '.nc', 
                             recursive = T, full.names = T)
    
    # mask grids to SIF vectors
    mask.file  <- tsif.files[basename(tsif.files) == 'SIFgrid.nc']
    tsif.files <- tsif.files[basename(tsif.files) != 'SIFgrid.nc']

    cat('grab.tsif(): Read mask file...\n')
    mask.dat <- nc_open(mask.file)
    #tsif.grd <- ncvar_get(mask.dat, 'mask')
    mask.lat <- ncvar_get(mask.dat, 'lat')
    mask.lon <- ncvar_get(mask.dat, 'lon')
    mask.id  <- ncvar_get(mask.dat, 'id')
    mask.grd <- expand.grid(lon = as.numeric(mask.lon), lat = as.numeric(mask.lat))
    mask.grd$id <- 1:nrow(mask.grd)
    sel.grd <- mask.grd %>% filter(id >= min(mask.id), id <= max(mask.id))

    # select TROPOMI SIF files
    tsif.timestr <- paste0(strsplit.to.df(gsub('.nc', '', basename(tsif.files)))$V3)
    if (!is.null(nhrs)) {  # if multiple nc files 
      time.info <- get.doy.from.timestr(paste0(timestr, '00'), nday = 1, nhrs = nhrs)
      tsif.file <- tsif.files[tsif.timestr >= substr(time.info$min.timestr, 1, 8) & 
                              tsif.timestr <= substr(time.info$max.timestr, 1, 8)]
    } else tsif.file <- tsif.files[tsif.timestr == timestr]

    # Start reading files
    if (length(tsif.file) == 0) {
        cat(paste('grab.tsif(): No downscaled TROPOMI SIF file found for', timestr, 
                  ', please check...\n'))
        return(NULL)

    } else {
      
      # form SIF into data.frame
      cat(paste('grab.tsif(): Need to process', length(tsif.file), 
                'file(s)...it takes a while\n'))
      
      rds.file <- file.path(dirname(tsif.file[1]), 'TSIF_JJA_mean.rds')
      if (length(tsif.file) == 92 & file.exists(rds.file)) {
        tsif.mean.df <- readRDS(rds.file) 

      } else {
        
        # calculate the sum
        tsif.sum <- list(ncvar_get(nc_open(tsif.file[1]), 'SIF'))

        if (length(tsif.file) > 1) {
          for (t in 2 : length(tsif.file)) {
            print(t)
            tsif.tmp <- list(ncvar_get(nc_open(tsif.file[t]), 'SIF'))
            tsif.sum <- mapply('+', tsif.sum, tsif.tmp, SIMPLIFY = FALSE)
          } # end for t

          # calculate the mean
          cat('grab.tsif: Cal-ing the mean SIF from multiple files...\n')
          tsif.mean <- as.numeric(unlist(tsif.sum)) / length(tsif.file)

        } else tsif.mean = tsif.sum # end if

        # add lat/lon id
        tsif.mean.df <- data.frame(SIF = tsif.mean, id = as.numeric(mask.id))
        if (length(tsif.file) == 92) saveRDS(tsif.mean.df, file = rds.file)
      } # end if


      # for debugging
      # minlon = -125; maxlon = -120; minlat = 35; maxlat = 40
      # minlon = -112.5; maxlon = -110.5; minlat = 39.5; maxlat = 41.5

      # select spatial regions
      cat('grab.tsif: forming TROPOMI SIF from vectors to grids...\n')
      crop.grd <- sel.grd %>% filter(lon >= minlon, lon <= maxlon, 
                                     lat >= minlat, lat <= maxlat)
      
      crop.tsif <- tsif.mean.df %>% filter(id >= min(crop.grd$id), 
                                           id <= max(crop.grd$id)) %>%
                                    left_join(crop.grd, by = 'id') %>% na.omit() 

      if (form == 'raster') {
        cat('grab.tsif: forming TROPOMI SIF from grids to rasterLayer...\n')
        crop.tsif <- suppressMessages(df2raster(crop.tsif[, c('SIF', 'lon', 'lat')]))
        #levelplot(tsif.rt, at = seq(0, 3, 0.1), main = timestr)
      } # end if
      
      return(crop.tsif)
    } # end if 

} # end of subroutine
