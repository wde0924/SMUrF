#' script to transform daily mean Reco to the correct form for predicting NEE
#' @author: Dien Wu, 08/09/2019 

#' input:
#' @param timestr YYYYMMDDHH
#' @param reco.file nc file containing rasterStack form of daily mean reco 
#'                 for each lat lon (see main_script_reco)
#' @param nhrs numbers of hours for aggregating/avergaing reco 
#' @param site.ext extent(minlon, maxlon, minlat, maxlat) with raster package

#' @return return mean reco within during @param nhrs since @param timestr


crop.smurf.reco <- function(timestr, reco.file, nhrs, site.ext = NULL) {

    reco.mean <- stack(reco.file, varname = 'Reco_mean')
    reco.sd <- stack(reco.file, varname = 'Reco_sd')

    if (!is.null(site.ext)) {
        reco.mean <- raster::crop(reco.mean, site.ext)
        reco.sd <- raster::crop(reco.sd, site.ext)    
    }

    # get reco rasterStack dimensions, second since 1970-01-01
    reco.names <- as.numeric(gsub('X', '', names(reco.mean)))
    reco.dates <- as.POSIXct(reco.names, origin = '1970-01-01 00:00:00', tz = 'UTC')
    reco.timestr <- as.numeric(paste0(format(reco.dates, '%Y%m%d')))
    reco.nd <- reco.timestr[2] - reco.timestr[1]       # temporal res in days


    ## select reco that falls within 'nhrs' from 'timestr'
    # nday for day interval, e.g., the temporal res of reco estiamtes, every 4 days
    if (nchar(timestr) != 10) stop('crop.smurf.reco(): "timestr" not in form of YYYYMMDDHH...\n')
    reco.doy.df  <- get.doy.from.timestr(timestr, nday = reco.nd, nhrs)
    min.timestr <- as.numeric(substr(reco.doy.df$min.timestr, 1, 8))
    max.timestr <- as.numeric(substr(reco.doy.df$max.timestr, 1, 8))
    layer.indx  <- which(substr(reco.timestr, 1, 8) >= min.timestr & 
                         substr(reco.timestr, 1, 8) <= max.timestr)
    
    # spatially crop and temporally select the rasterStack
    reco.mean.rt <- raster::subset(reco.mean, layer.indx)
    reco.sd.rt <- raster::subset(reco.sd, layer.indx)


    ## calculate mean reco if there are several 4 or 8 days layers for backward duration 
    # in other words, length(layer.indx) > 1 
    if (length(layer.indx) > 1) {
        reco.mean.rt <- mean(reco.mean.rt); reco.sd.rt <- mean(reco.sd.rt)
    }

    # combine reco mean and SD into one rasterStack 
    reco.stk <- stack(reco.mean.rt, reco.sd.rt)
    names(reco.stk) <- c('Reco_mean', 'Reco_sd')

    return(reco.stk)
}