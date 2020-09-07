#' script to preprocess Tair from Daymet, DW, 05/22/2019 
#' @param daymet.path path for daymet
#' @param daymet.files daymet filenames separate by years, e.g., daymet_v3_tmax_2010_na.nc4
#'               filename needs to contain parameters, e.g., tmax, tmin
#' @param timestr YYYYMMDDHH
#' @param ext need to be generated from extent() in raster package; 
#'            i.e., extent(minlon, maxlon, minlat, maxlat)
#' @param nhrs numbers of hours backward or forward for potential X-STILT runs, can be NULL
#' @param tz timezone of `timestr`, default is UTC

#' @return rasterLayer of daily mean Temp from daymet (calcualte from Tmin, Tmax)

# daymet.file only contains data for a single year !!
prep.tair.daymet <- function(tair.path, tair.varname = 'daymet_v3', timestr, ext, 
                             nhrs = NULL, tz = 'UTC') {

    library(raster)

    # load all daily max and min temp files 
    daymet.files <- list.files(tair.path, tair.varname, full.names = T, recursive = T)
    daymet.files <- daymet.files[grepl(substr(timestr, 1, 4), basename(daymet.files))]
    if (length(daymet.files) == 0) 
        stop(paste0('prep.tair.daymet(): NO daymet file containing "', 
                    substr(timestr, 1, 4), '" found\n'))

    # subset days, nday = 1 since tair is daily 
    tair.doy.df <- get.doy.from.timestr(timestr, nday = 1, nhrs = nhrs, tz = tz)
    tair.band   <- seq(tair.doy.df$min.doy, tair.doy.df$max.doy)

    # https://cran.r-project.org/web/packages/daymetr/vignettes/daymetr-vignette.html
    # read as RasterStack or RasterLayer
    if (FALSE %in% file.exists(daymet.files)) 
        stop('prep.tair.daymet(): NO tmax or tmin file, please check..\n')
    tmax.stk <- raster::stack(daymet.files[grep('tmax', daymet.files)], bands = tair.band)
    tmin.stk <- raster::stack(daymet.files[grep('tmin', daymet.files)], bands = tair.band)
    #tmax.stk <- raster::readAll(tmax.stk)
    cat(paste('prep.tair.daymet(): Done loading all tmax and tmin for North America',
              'and start to average Tmin and Tmax for daily mean T...\n'))

    # need to crop the data before calculating the Tmean, calculate the extent in meters
    daymet.ext   <- get.daymet.ext(ext)
    sel.tmax.stk <- raster::crop(tmax.stk, daymet.ext)
    sel.tmin.stk <- raster::crop(tmin.stk, daymet.ext)

    # average both tmin and tmax for ndays (as SIF is only available for 4 or 8 days)
    # and reproject to lat lon
    sel.tmean.rt <- mean(sel.tmin.stk, sel.tmax.stk)  

    # project to longlat coordinate with bilinear interpolation
    raster::projection(sel.tmean.rt) <- 
    '+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +a=6378137 +b=6356752.314706705 +units=m +no_defs'
    
    # reproject to longlat coordinates, with bilinear interpolation
    cat('prep.tair.daymet(): Reprojecting Tmean...\n')
    sel.tmean.pj <- raster::projectRaster(sel.tmean.rt, crs = '+init=epsg:4326')
    
    # further crop Tmean according to desired spatial region around site
    crop.tmean.pj <- raster::crop(sel.tmean.pj, ext)
    return(crop.tmean.pj)

} # end of subrotuine


# inner function to get daymet extent in meters from latlon coordinates
get.daymet.ext <- function(ext) {

    library(geosphere); library(raster)
    lcc.lat <- 42.5; lcc.lon <- -100

    # because there will be spatial mismatches between lcc and longlat projection
    # we allow for an additional 4x4 degree area
    # high latitude locations need more extention
    p2 <- data.frame(lon = c(ext[1] - 4, ext[2] + 4, lcc.lon, lcc.lon), 
                     lat = c(lcc.lat, lcc.lat, ext[3] - 4, ext[4] + 4))
    
    p2$dist <- distCosine(c(lcc.lon, lcc.lat), p2 = p2) # calculate dist in m

    # correct the sign of dist by comparing lat, lon
    p2 <- p2 %>% mutate(lcc.lon = lcc.lon, lcc.lat = lcc.lat, 
                        loc = lon + lat, lcc.loc = lcc.lon + lcc.lat, 
                        dist = ifelse(loc < lcc.loc, -dist, dist))

    extent(p2$dist)
} # end of subroutine











# get tile ID info using tile_outlines from daymetr package
get.tile.id <- function(ext) {

    library(daymetr); library(raster)
    tile.sp <- tile_outlines
    tile.df <- raster::as.data.frame(tile.sp) %>% group_by(TileID) %>% 
               mutate(lat = (YMin + Ymax) / 2, lon = (XMin + XMax) / 2) %>% ungroup()

    tile.rt <- df2raster(tile.df)
    sel.tile.rt <- crop(tile.rt, ext)
    sel.tile.df <- raster::as.data.frame(sel.tile.rt, xy = T)
    #tile.id <- sel.tile.df$TileID
    return(sel.tile.df)
}


