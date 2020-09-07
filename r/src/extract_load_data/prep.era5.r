#' script to preprocess hourly Tair from ERA5, DW, 06/27/2019 
#' return hourly Tair as rasterBrick (if multiple hours) or rasterLayer (if one single hour)
#' no daily mean calculation performed yet

#' ERA5 data is downloaded from
#' https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
#' due to large file size, store 2m Tair (K) for each month separately 
#' *** store filename as 'AT2m_YYYYMM.nc' e.g., AT2m_201401.nc 

#' UTC reported in ERA5 represents the end hour!!!, DW, 05/26/2020
#' see https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation
#' so I need to move the UTC from ENDING hour to STARTING hour 
#     to be consistent with the hour convention for fluxes

#' @param timestr in the form of YYYYMMDDHH
#' @param tz indicating the time zone of @param timestr
#' @param era5.varname like '2T_', 'SSRD_', 'STL1_', 'SWVL1_' indicating the variables 
#'                     defined by ERA5 

prep.era5 <- function(era5.path, era5.varname = c('2T', 'SSRD', 'STL1', 'SWVL1')[1], 
                      timestr = '2018010100', ext = NULL, nhrs = NULL, tz = 'UTC') {
    
    library(raster); library(dplyr); library(ncdf4)
    if (nchar(timestr) != 10) stop('prep.era5(): timestr should be in form of YYYYMMDDHH')

    # load all Tair files
    era5.files <- list.files(era5.path, era5.varname, full.names = T, recursive = T)

    # ------------------------------------------------------------------------ #
    # locate the files needed for `timestr` and `nhrs` (if nhrs is not NULL)
    doy.df <- get.doy.from.timestr(timestr, nday = 1, nhrs = nhrs, tz = tz)
    seq.dates <- seq(as.POSIXct(doy.df$min.timestr, tz = tz, format = '%Y%m%d%H'), 
                     as.POSIXct(doy.df$max.timestr, tz = tz, format = '%Y%m%d%H'), 
                     by = 1 * 60 * 60)    # by in second
    if (!is.null(nhrs)) seq.dates <- seq.dates[1 : nhrs]

    # if the last element in "seq.dates" is the last day of a month, 
    # need to grab an extra nc file, DW, 05/26/2020 
    last.date <- seq.dates[length(seq.dates)]

    # check to see if last.date is the last day of a month
    dom <- lubridate::days_in_month(as.numeric(substr(last.date, 6, 7)))
    if (as.numeric(substr(last.date, 9, 10)) == dom && substr(last.date, 12, 13) == '23') {
        cat('prep.era5(): encountered the last day of a month, need to grab an additional file\n')
        new.date <- last.date + 1 * 3600 
        add.yyyymm <- format(new.date, format = '%Y%m')
    } else add.yyyymm <- NULL 

    uni.yyyymm <- c(unique(format(seq.dates, format = '%Y%m')), add.yyyymm)


    # search for monthly ERA5 files first (hourly Tair stored in a montly manner)
    # if not, look for annual ERA5 file (hourly Tair stored in a yearly manner)
    era5.file <- NULL
    for (m in uni.yyyymm) {
        tmp.file <- era5.files[grepl(m, era5.files)]

        if (length(tmp.file) == 0) {   # if no monthly file, look for annual file
            tmp.file <- era5.files[grepl(substr(timestr, 1, 4), era5.files)]
            if (length(era5.file) == 0) 
                stop(paste('prep.era5(): NO era5 file found for', m))
        }   # end if
        
        era5.file <- c(era5.file, tmp.file)
    }


    # ------------------------------------------------------------------------ #
    # read all required hourly Tair, convert deg K to deg C
    cat(paste('prep.era5(): Reading hourly', era5.varname, '...\n'))

    era5.stk <- NULL
    for (f in 1 : length(era5.file)) {
       
        tmp.times <- ncvar_get(nc_open(era5.file[f]), 'time')
        tmp.end.dates <- as.POSIXct(tmp.times * 3600, origin = '1900-01-01 00:00:00', 
                                    tz = 'UTC')

        # !!!! hours in ERA5 stands for the ending hours
        # so we now move it to the starting hour, DW, 05/26/2020 
        tmp.start.dates <- tmp.end.dates - 1 * 3600

        # locate the correct time period, read Tair as rasterStack
        tmp.bands <- which(tmp.start.dates %in% seq.dates)
        tmp.stk <- raster::stack(era5.file[f], bands = tmp.bands)   

        # the dates that raster processed are WRONG!!!
        names(tmp.stk) <- tmp.start.dates[tmp.bands] 

        # store data 
        if (f == 1) era5.stk <- tmp.stk
        if (f >  1) era5.stk <- stack(era5.stk, tmp.stk)
    }   # end for f


    # ------------------------------------------------------------------------ #
    # adjust ERA5's longitude from 0-360 to -180 to 180 
    # ERA5 reports temperatures right at each lat/lon, not for a tile, 
    if (extent(era5.stk)[2] > 181) {
        extent(era5.stk) <- extent(0, 360, -90.125, 90.125)
        east <- crop(era5.stk, extent(0, 180, -90.125, 90.125))
        west <- crop(era5.stk, extent(180, 360, -90.125, 90.125))

        # then change extent of west to negative long
        extent(west) <- c(-180, 0, -90.125, 90.125)
        era5.ref.stk <- merge(west, east)
        names(era5.ref.stk) <- names(era5.stk)
    } else { era5.ref.stk <- era5.stk }

    # crop the data according to spatial ext (can be generated from extent())
    if (!is.null(ext)) era5.ref.stk <- crop(era5.ref.stk, ext)
    
    # ------------------------------------------------------------------------ #
    # unit corrections for '2T', 'SSRD', 'STL1', 'SWVL1'    
    if (grepl('2T', era5.varname) | grepl('STL', era5.varname)) {

        cat('prep.era5(): *** unit conversion (degC for temps)\n')
        era5.ref.stk <- calc(era5.ref.stk, fun = function(x) { x - 273.15 }) 

    } else if (era5.varname == 'SSRD'){

        # if ssrd is negative like -4.656613e-10, force it to be zero 
        # also convert J/m2 to W/m2, dividing hrly total radiation by 3600 s
        cat('prep.era5(): *** unit conversion (W/m2 for radiation)\n')
        era5.ref.stk[era5.ref.stk < 0] <- 0
        era5.ref.stk <- era5.ref.stk / 3600 # now in W/m2

    } else era5.ref.stk <- era5.ref.stk

    # rasterVis::levelplot(era5.ref.stk)
    return(era5.ref.stk)
}   # end of subroutine