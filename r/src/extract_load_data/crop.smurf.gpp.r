# script to transform daily mean GPP to the correct form for predicting Reco 
#' @author: Dien Wu, 05/28/2019 

#' input:
#' @param timestr YYYYMMDDHH
#' @param gpp.file nc file containing rasterStack form of daily mean GPP 
#'                 for each lat lon (see main_script_GPP)
#' @param nhrs numbers of hours for aggregating/avergaing GPP 
#' @param site.ext extent(minlon, maxlon, minlat, maxlat) with raster package

#' @return return mean GPP within during @param nhrs since @param timestr

#' updates:
#' 08/01/2019, read daily GPP from nc file

crop.smurf.gpp <- function(timestr, gpp.file, nhrs, site.ext = NULL, 
                           time.integrate = T) {

    gpp.mean <- stack(gpp.file, varname = 'GPP_mean')
    gpp.sd <- stack(gpp.file, varname = 'GPP_sd')

    if (!is.null(site.ext)) {
        gpp.mean <- raster::crop(gpp.mean, site.ext)
        gpp.sd <- raster::crop(gpp.sd, site.ext)    
    }

    # get GPP rasterStack dimensions, second since 1970-01-01
    gpp.names <- as.numeric(gsub('X', '', names(gpp.mean)))
    gpp.dates <- as.POSIXct(gpp.names, origin = '1970-01-01 00:00:00', tz = 'UTC')
    gpp.timestr <- as.numeric(paste0(format(gpp.dates, '%Y%m%d')))
    gpp.nd <- gpp.timestr[2] - gpp.timestr[1]       # temporal res in days

    ## select GPP that falls within 'nhrs' from 'timestr'
    # nday for day interval, e.g., the temporal res of GPP estiamtes, every 4 days
    if (nchar(timestr) != 10) stop('crop.smurf.gpp(): "timestr" not in form of YYYYMMDDHH...\n')
    gpp.doy.df  <- get.doy.from.timestr(timestr, nday = gpp.nd, nhrs)
    min.timestr <- as.numeric(substr(gpp.doy.df$min.timestr, 1, 8))
    max.timestr <- as.numeric(substr(gpp.doy.df$max.timestr, 1, 8))
    layer.indx  <- which(substr(gpp.timestr, 1, 8) >= min.timestr & 
                         substr(gpp.timestr, 1, 8) <= max.timestr)
    cat(paste('crop.smurf.gpp(): will load files (every', gpp.nd, 
              'day(s)) with timestr from', min.timestr, 'to', max.timestr, 
              ';total', length(layer.indx), 'of rasterLayers\n'))

    # spatially crop and temporally select the rasterStack
    gpp.mean.rt <- raster::subset(gpp.mean, layer.indx)
    gpp.sd.rt <- raster::subset(gpp.sd, layer.indx)

    ## calculate mean GPP if there are several 4 or 8 days layers for backward duration 
    # in other words, length(layer.indx) > 1 
    if (length(layer.indx) > 1 & time.integrate == T) {
        gpp.mean.rt <- mean(gpp.mean.rt); gpp.sd.rt <- mean(gpp.sd.rt)
        
        # combine GPP mean and SD into one rasterStack 
        gpp.stk <- stack(gpp.mean.rt, gpp.sd.rt)
        names(gpp.stk) <- c('GPP_mean', 'GPP_sd')
        return(gpp.stk)

    } else if ( length(layer.indx) > 1 & time.integrate == F) {
        gpp.list <- list(GPP_mean = gpp.mean.rt, GPP_sd = gpp.sd.rt)
        return(gpp.list)

    } else {    # only one layer
        
        # combine GPP mean and SD into one rasterStack 
        gpp.stk <- stack(gpp.mean.rt, gpp.sd.rt)
        names(gpp.stk) <- c('GPP_mean', 'GPP_sd')
        return(gpp.stk)
    }   # end if


}
