#' scripts to save several rasterStack into nc files

# --------------------------------------------------------------------------- #
#' @author Dien Wu
#' @param varnames a vector of the names of all variables, choose from 
#'                 'SIF_mean', 'SIF_sd', 'GPP_mean', 'GPP_sd', 'Reco_mean', 
#'                 'Reco_sd', 'NEE_mean', 'NEE_sd'
#' @param varunits a vector of units of all variables, should match @param varnames
#' @param longnames a vector of the long descriptions of all variables
#'                  should match @param varnames
#' @param stk.list list of rasterStacks whose list names 
#'                 should match the variables in @param varnames, 

#'        ***** names of each rasterLayer or stack should indicate the time, 
#'              with correct form in POSIXct() *****

#'        ***** ALL rasterStacks in this list should have the SAME extent, 
#'                  lat/lon/time dimensions ***** 

#' @param timescale whether daily mean or hourly mean fluxes
#'                  will reflect in the longname of the variable in nc file
#' @param zformat the format of z values (time indicators), e.g., 
#'                'X%Y.%m.%d' (year, month, day) OR 'X%Y.%m.%d.%H.%M.%S' (accuracy to second)
#'        ***** if there is no time dimension in @param stk.list, set zformat = NA

save.raster2nc <- function(varnames = c('SIF_mean', 'SIF_sd', 'GPP_mean', 'GPP_sd', 
                                        'Reco_mean', 'Reco_sd', 'NEE_mean', 'NEE_sd'), 
                           varunits = c(rep('mW m−2 nm−1 sr−1', 2), rep('umol m-2 s-1', 6)), 
                           longnames, 
                           zformat = c('X%Y.%m.%d', 'X%Y.%m.%d.%H.%M.%S', NA)[1], 
                           stk.list = list(), filename) {

    library(raster); library(ncdf4)
    if ( is.null(names(stk.list))) stop('save.raster2nc: need to assign names to stk.list')
    if ( TRUE %in% (names(stk.list) != varnames) ) 
        stop('save.raster2nc: need to assign names to stk.list that match varnames')

    # construct a data frame 
    nc.df <- data.frame(varname = varnames, varunit = varunits, 
                        longname = longnames, zformat = zformat, 
                        stringsAsFactors = F)

    # initialize with False, for storing the first var, 
    # overwrite and create new nc file, rather than append for adding var
    append <- F; nc <- NULL 
    for ( n in 1 : nrow(nc.df) ) {

        if (n > 1) append <- T  # add variables if not the 1st round
        tmp.df <- nc.df[n, ]
        tmp.stk <- stk.list[[ which(names(stk.list) == tmp.df$varname) ]]

        nc <- ncdim.ncput(stk = tmp.stk, varname = tmp.df$varname, 
                          varunit = tmp.df$varunit, longname = tmp.df$longname, 
                          zformat = tmp.df$zformat, append, filename, nc)

    }   # end for n

    if (!is.na(zformat)) {
        ncatt_put(nc, 0, 'crs', '+proj=longlat')
        ncatt_put(nc, 0, 'crs_format', 'PROJ.4')
    }

    ncatt_put(nc, 0, 'documentation', 'github.com/uataq/SMUrF')
    ncatt_put(nc, 0, 'title', 'SMUrF Biospheric Fluxes')
    ncatt_put(nc, 0, 'time_created', format(Sys.time(), tz = 'UTC'))
    nc_close(nc)

    cat(paste('save.raster2nc(): nc created in', filename, '\n'))
}   # end of script


# --------------------------------------------------------------------------- #
# need an inner function ------------------------------------------------- #
# parts of the code are inspired from 'write_footprint.r', written by Ben Fasoli
# edited by Dien Wu, 07/25/2019 

# stk could be rasterLayer or rasterStack for mean and SD of the fluxes
# if append == T, nc needs to have values, will skip the portion that stores lat/lon
# if append == F, first time writing a file and lat/lon 

# 08/05/2019 bug fixed for converting rasterLayer to 2D matrix, 
#            transpose in raster::as.array won't work, DW
# zformat for the format of z values (time indicators)

ncdim.ncput <- function(stk, varname = 'GPP', varunit = 'umol m-2 s-1', 
                        longname = 'Daily mean Gross Primary Production', 
                        zformat = c('X%Y.%m.%d', 'X%Y.%m.%d.%H.%M.%S', NA)[1], 
                        append = F, filename = NULL, nc = NULL) {

    # determine lat/lon info from rasterStack
    xres <- xres(stk); yres <- yres(stk)
    xmn <- extent(stk)[1]; xmx <- extent(stk)[2]
    ymn <- extent(stk)[3]; ymx <- extent(stk)[4]
    
    # set grid of fluxes, increasing trend, these are lower left lat/lon
    glong <- head(seq(xmn, xmx, by = xres), -1)
    glati <- head(seq(ymn, ymx, by = yres), -1)

    # xy dimensions in lat/lon or alternative projection. 
    # move lower left to centered lat/lon
    xdim <- ncdim_def('lon', 'degrees_east',  glong + xres/2)
    ydim <- ncdim_def('lat', 'degrees_north', glati + yres/2)

    # -999 for missing
    fvar <- ncvar_def(varname, varunit, list(xdim, ydim), -999) 

    # if names of each rasterLayer or stack indicate the time, 
    # with correct form in POSIXct() initially
    if (!is.na(zformat)) {
        options(scipen = 999)
        gname <- names(stk)     # names of rasterLayer/Stack, indicating time
        time_out <- as.POSIXct(gname, format = zformat, tz = 'UTC')
        tdim <- ncdim_def('time', 'seconds since 1970-01-01 00:00:00Z',
                                   as.numeric(time_out))
        # -999 for missing
        fvar <- ncvar_def(varname, varunit, list(xdim, ydim, tdim), -999) 
    }

    # convert raster stack to array and -------------------------------------- #
    # flip latitude from descending to ascending trend (raster uses descending)

    # *** raster::as.array transpose only works for 3D array
    # use t() for 2D matrix, DW, 08/05/2019
    if (nlayers(stk) == 1) {    # dim of raster: [lat, lon...]

        fmatrix <- t(raster::as.matrix(stk) )  # [lon, lat]
        fflip <- fmatrix[, length(glati):1]   # flip lat, now in ascending order

    } else {
        
        # rasterStack
        farray <- raster::as.array(stk, transpose = T)   # [lon, lat, time]
        fflip <- farray[, length(glati):1, ]  # flip lat, now in ascending order
    }   # end if

    # if append == T, nc needs to have values, will skip the portion that stores lat/lon
    # if append == F, first time writing a file and lat/lon 
    if (append == F) {      # Projection specific xy definitions
        nc <- nc_create(filename, list(fvar), force_v4 = append)
        ncatt_put(nc, 'lon', 'standard_name', 'longitude')
        ncatt_put(nc, 'lon', 'long_name', 'longitude at cell center')
        ncatt_put(nc, 'lat', 'standard_name', 'latitude')
        ncatt_put(nc, 'lat', 'long_name', 'latitude at cell center')

        # Insert flux data
        ncvar_put(nc, fvar, fflip)
        ncatt_put(nc, varname, 'standard_name', tolower(varname))
        ncatt_put(nc, varname, 'long_name', longname)

    } else {
        
        # add another set of fluxes to existing files
        nc <- ncvar_add(nc, fvar) 
        ncvar_put(nc, fvar, fflip)
        ncatt_put(nc, varname, 'standard_name', tolower(varname))
        ncatt_put(nc, varname, 'long_name', longname)
    }  # end if append
   
    if (!is.na(zformat)) {
        ncatt_put(nc, 'time', 'standard_name', 'time')
        ncatt_put(nc, 'time', 'long_name', 'utc time')
        ncatt_put(nc, 'time', 'calendar', 'standard')
    }

    return(nc)
}   # end of ncdim.ncput()
