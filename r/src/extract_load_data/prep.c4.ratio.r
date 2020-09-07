# function to prepare C3-C4 ratio and disaggregate to MODIS land cover 
# Dien Wu, 03/29/2020 

prep.c4.ratio <- function(smurf_wd, lc.path, lc.pattern, yr, lc.max.yr, reg.name, 
                          minlon, maxlon, minlat, maxlat, gpp.path) {
    
    out.fn <- file.path(gpp.path, paste0('C4_ratio_', reg.name, '.tif'))

    if (file.exists(out.fn)) {
        return(out.fn)
        
    } else {
        
        reg.ext <- extent(minlon, maxlon, minlat, maxlat)

        # prepare MCD12 IGBP land cover
        lc.rt <- prep.mcd12(lc.path, lc.pattern, yr, lc.max.yr, reg.name, reg.ext)

        ratio.fn <- file.path(smurf_wd, 'data/C4_relative_fraction.tif')
        if (!file.exists(ratio.fn)) 
            stop('NO nc file found for C3-C4 ratio..please check repository\n')
        c4.rt <- crop(raster(ratio.fn, varname = 'C4_rel_frac'), extent(lc.rt))

        # project to 500 m 
        fact <- round(res(c4.rt)[1] / res(lc.rt)[1])
        c4.ds <- projectRaster(raster::disaggregate(c4.rt, fact = fact), lc.rt)
        
        writeRaster(c4.ds, out.fn, overwrite = T, format = 'GTiff', 
                    varname = 'C4_rel_frac', xname = 'lon', yname = 'lat', 
                    longname = 'reprojected 500 m C4 relative fractions')

        return(out.fn)

    }   # end if

}   # end of function