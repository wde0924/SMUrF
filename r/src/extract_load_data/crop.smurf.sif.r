# script to transform daily mean sif to the correct form for predicting Reco 
#' @author: Dien Wu, 05/28/2019 

#' input:
#' @param timestr YYYYMMDDHH
#' @param sif.file nc file containing rasterStack form of daily mean sif 
#'                 for each lat lon (see main_script_sif)
#' @param nhrs numbers of hours for aggregating/avergaing sif 
#' @param site.ext extent(minlon, maxlon, minlat, maxlat) with raster package

#' @return return mean sif within during @param nhrs since @param timestr

#' updates:
#' 08/01/2019, read daily sif from nc file

crop.smurf.sif <- function(timestr, sif.file, nhrs, site.ext = NULL) {

    sif.mean <- stack(sif.file, varname = 'SIF_mean')
    if (!is.null(site.ext)) sif.mean <- raster::crop(sif.mean, site.ext)

    # get sif rasterStack dimensions, second since 1970-01-01
    sif.names <- as.numeric(gsub('X', '', names(sif.mean)))
    sif.dates <- as.POSIXct(sif.names, origin = '1970-01-01 00:00:00', tz = 'UTC')
    sif.timestr <- as.numeric(paste0(format(sif.dates, '%Y%m%d')))
    sif.nd <- sif.timestr[2] - sif.timestr[1]       # temporal res in days


    ## select sif that falls within 'nhrs' from 'timestr'
    # nday for day interval, e.g., the temporal res of sif estiamtes, every 4 days
    if (nchar(timestr) != 10) stop('crop.smurf.sif(): "timestr" not in form of YYYYMMDDHH...\n')
    sif.doy.df  <- get.doy.from.timestr(timestr, nday = sif.nd, nhrs)
    min.timestr <- as.numeric(substr(sif.doy.df$min.timestr, 1, 8))
    max.timestr <- as.numeric(substr(sif.doy.df$max.timestr, 1, 8))
    layer.indx  <- which(substr(sif.timestr, 1, 8) >= min.timestr & 
                         substr(sif.timestr, 1, 8) <= max.timestr)
    cat(paste('crop.smurf.sif(): will load files (every', sif.nd, 
              'day(s)) with timestr from', min.timestr, 'to', max.timestr, 
              ';total', length(layer.indx), 'of rasterLayers\n'))

    # spatially crop and temporally select the rasterStack
    sif.mean.rt <- raster::subset(sif.mean, layer.indx)

    ## calculate mean sif if there are several 4 or 8 days layers for backward duration 
    # in other words, length(layer.indx) > 1 
    if (length(layer.indx) > 1) sif.mean.rt <- mean(sif.mean.rt)

    # combine sif mean and SD into one rasterStack 
    names(sif.mean.rt) <- 'SIF_mean'
    return(sif.mean.rt)
}
