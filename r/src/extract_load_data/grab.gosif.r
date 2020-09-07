# subroutine to grab GOSIF from Li and Xiao, 2019, MDPI
# by Dien Wu, 03/21/2019

# if year > 2017, use GOSIF for year 2017
# timestr: in the form of YYYYMMDDHH in UTC, HH for overpass hour
grab.gosif <- function(gosif.path, timestr = 2015010100, lon.lat = NULL, 
                       ext = NULL, nhrs = NULL, 
                       output.form = c('raster', 'df')[1]) {

    # call a funcion to get DOY from YYYYMMDD, given day intervals
    time.info <- get.doy.from.timestr(timestr, nday = 8, nhrs = nhrs)
    search.timestr <- time.info$find.yrdoy

    # if year > 2017, use gosif for year 2017, DW, 01/10/2019
    if (time.info$yr0 > 2017) {
      cat('grab.gosif(): *** request year > 2017, use data of year 2017 instead...\n')
      search.timestr <- paste0('2017', substr(search.timestr, 5, 10))
    } 

    # select the ones that need for scaling GPP 
    gosif.files <- list.files(path = gosif.path, pattern = '.tif', 
                              recursive = T, full.names = T)
    gosif.file <- gosif.files[grep(search.timestr, gosif.files)]

    if (length(gosif.file) == 0) {
        cat(paste('grab.gosif(): No gosif file found for', search.timestr, 
                  ', please check...\n'))
        return()

    } else {

        # grab instantaneous and daily averaged SIF
        gosif <- raster(gosif.file)
        extent(gosif) <- extent(-180, 180, -90, 90)

        # 100 for NA, assign as NA for water
        gosif[gosif == 100] <- NA
    
        # crop spatial gosif for 20x20deg footprint domain
        if (!is.null(lon.lat) & is.null(ext)) ext <- extent(lon.lat$minlon, 
                                                            lon.lat$maxlon, 
                                                            lon.lat$minlat, 
                                                            lon.lat$maxlat)
        sel.gosif <- crop(gosif, ext)
        
        if (output.form == 'df') {
            sel.gosif <- raster::as.data.frame(sel.gosif, xy = T)
            colnames(sel.gosif) <- list('lon', 'lat', 'GOSIF')
        }
        
        return(sel.gosif)
    } # end if 

} # end of subroutine
