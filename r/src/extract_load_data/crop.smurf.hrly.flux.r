#' script to transform daily mean Reco to the correct form for predicting flux
#' @author: Dien Wu, 08/09/2019 

#' input:
#' @param timestr YYYYMMDDHH
#' @param flux.file nc file containing rasterStack form of daily mean reco 
#'                 for each lat lon (see main_script_reco)
#' @param nhrs numbers of hours for aggregating/avergaing reco 
#' @param site.ext extent(minlon, maxlon, minlat, maxlat) with raster package
#' @param hr if !is.null(hr) == T, only grab the flux during this hour interval
#'           e.g., hr = 4, grab hourly mean fluxes during 04UTC to 05UTC
#' @return return mean reco within during @param nhrs 
#'         (including the one hourly flux starting at timestr) from @param timestr


# if intTF == T, return mean of hourly fluxes, if not, return all hourly fluxes
crop.smurf.hrly.flux <- function(timestr, flux.file, site.ext = NULL,
                                 varname = c('GPP_mean', 'Reco_mean', 'NEE_mean')[3], 
                                 nhrs = NULL, hr = NULL, intTF = F, 
                                 res = c(1/120, 1/20)[2]) {

    if ( res == 1/120 ) {   # high-res SMUrF, flux at 1km

        doy.df <- get.doy.from.timestr(timestr, nday = 1 / 24, nhrs = nhrs)
        file.timestr <- strsplit.to.df(gsub('.nc', '', basename(flux.file)))$V4
        sel.flux.file <- flux.file[file.timestr >= substr(doy.df$min.timestr, 1, 8) & 
                                   file.timestr <= substr(doy.df$max.timestr, 1, 8)]
        
        flux.mean <- stack(sel.flux.file[1], varname = varname)
        if (!is.null(site.ext)) flux.mean <- raster::crop(flux.mean, site.ext)
        if (length(sel.flux.file) > 1) {
            for (f in 2 : length(sel.flux.file)) { 
                print(f)
                tmp.flux <- stack(sel.flux.file[f], varname = varname)
                if (!is.null(site.ext)) tmp.flux <- raster::crop(tmp.flux, site.ext)
                flux.mean <- stack(flux.mean, tmp.flux)
            }   # end for f
        }   # end if

    } else {
        flux.mean <- stack(flux.file, varname = varname)
        if (!is.null(site.ext)) flux.mean <- raster::crop(flux.mean, site.ext)
    }

    # string must contain 'e' for scientific notation
    fix.sci.not <- function(string) as.character(as.integer(gsub('e.', 'e+', string)))

    # get reco rasterStack dimensions, second since 1970-01-01
    flux.names <- gsub('X', '', names(flux.mean))
    flux.names[grepl('e.', flux.names)] <- fix.sci.not(flux.names[grepl('e.', flux.names)])
    flux.names <- as.numeric(flux.names)
    flux.dates <- as.POSIXct(flux.names, origin = '1970-01-01 00:00:00', tz = 'UTC')
    flux.timestr <- as.numeric(paste0(format(flux.dates, '%Y%m%d%H')))
    flux.nhr <- flux.timestr[2] - flux.timestr[1]       # temporal res in days

    ## select reco that falls within 'nhrs' from 'timestr'
    # nday for day interval, e.g., the temporal res of reco estiamtes, every 4 days
    if (nchar(timestr) != 10) 
       stop('crop.smurf.hrly.flux(): "timestr" not in form of YYYYMMDDHH...\n')
    
    flux.doy.df <- get.doy.from.timestr(timestr, nday = flux.nhr / 24, nhrs = nhrs)
    min.timestr <- as.numeric(flux.doy.df$min.timestr)
    max.timestr <- as.numeric(flux.doy.df$max.timestr)
    layer.indx  <- which(flux.timestr >= min.timestr & flux.timestr <= max.timestr)

    if (!is.null(hr)) layer.indx <- which(substr(flux.timestr, 1, 8) >= min.timestr & 
                                          substr(flux.timestr, 1, 8) <= max.timestr & 
                                          substr(flux.timestr, 9, 10) == hr)
    
    # spatially crop and temporally select the rasterStack
    flux.stk <- raster::subset(flux.mean, layer.indx)
    if (nlayers(flux.stk) > 1 && intTF) {
        flux.mean.rt <- mean(flux.stk); names(flux.mean.rt) <- varname
        return(flux.mean.rt)

    } else {
        names(flux.stk) <- flux.timestr[layer.indx]
        return(flux.stk)
    }   # end if

}   # end of subroutine