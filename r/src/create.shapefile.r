# subroutine to create shapefile around a city, given rectangle domain
# DW, 03/18/2019

# add site.ext: extent(minlon, maxlon, minlat, maxlat)
create.shapefile <- function(site, site.ext = NULL, dlon = 20, dlat = 20, 
                             outpath, overwrite = F) {

    if (is.null(site.ext)) {
        lon.lat <- get.lon.lat(site, dlon, dlat) 
        coor <- data.frame(
                lon = c(rep(lon.lat$minlon, 2), rep(lon.lat$maxlon, 2)), 
                lat = c(lon.lat$minlat, rep(lon.lat$maxlat, 2), lon.lat$minlat)
                )
    } else {
        coor <- data.frame(lon = c(rep(site.ext[1], 2), rep(site.ext[2], 2)), 
                           lat = c(site.ext[3], rep(site.ext[4], 2), site.ext[3]))

    }  # end if
   
    tmp <- coor
    sr  <- Polygon(tmp)  # convert data frame to polygon class
    srs <- Polygons(list(sr), ID = 1)   # convert to polygonS class

    # create SpatialPolygons object
    sp <- SpatialPolygons(list(srs), proj4string = CRS('+proj=longlat +datum=WGS84'))
    shp.filename <- file.path(outpath, paste0(site, '_shapefile.shp'))

    # save as shapefiles
    if (!file.exists(shp.filename) | overwrite) 
        raster::shapefile(sp, shp.filename, overwrite = overwrite) 

    # create zip files
    files <- list.files(outpath, paste0(site, '_shapefile'), full.names = T)
    zipfile <- file.path(outpath, paste0(site, '_shapefile.zip'))
    zip(zipfile = zipfile, files)

    return(zipfile)
}
