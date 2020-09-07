# read above-ground biomass data, by Dien Wu, 04/17/2019
# GlobBiomass: above ground biomass (AGB, unit: tons/ha i.e., Mg/ha) for the year 2010
# Definition: the mass, expressed as oven-dry weight of the woody parts 
# (stem, bark, branches and twigs) of all living trees excluding stump and roots.

# ORNL: 1x1deg, several raster layers, need to insert 'agb.var', 
#       e.g., AGB_ha, AboveGroundBiomass, BrBCEF, ConBCEF, ConiferShare...
# ext: extent of raster layer, extent(xmin, xmax, ymin, ymax)
# whether get raster or data.frame form

# *** will store the reformed tif file and return the tif filename
grab.agb.GlobBiomass <- function(agb.path, minlon, maxlon, minlat, maxlat){

    # reform AGB 
    ext.str <- paste(minlon, maxlon, minlat, maxlat, sep = '_')
    out.fn <- file.path(agb.path, paste0('reform_AGB_', ext.str, '.tif'))

    if (!file.exists(out.fn)) {
    
        # find the correct 40x40deg tiles from globBiomass excluding error files
        # find the correct biomass file based on ext
        ext <- extent(minlon, maxlon, minlat, maxlat)
        lon.int <- seq(-180, 180, 40)
        lat.int <- seq(-40, 80, 40)
        lon.indx <- unique(c(findInterval(minlon, lon.int), 
                             findInterval(maxlon - 1E-10, lon.int)))
        lat.indx <- unique(c(findInterval(minlat, lat.int), 
                             findInterval(maxlat - 1E-10, lat.int))) + 1

        agb.lon <- lon.int[lon.indx]
        agb.lon.str <- ifelse(agb.lon < 0, 
                              paste0('W', formatC(abs(agb.lon), width = 3, flag = 0)), 
                              paste0('E', formatC(agb.lon, width = 3, flag = 0)))
        
        agb.lat <- lat.int[lat.indx]
        agb.lat.str <- ifelse(agb.lat < 0, 
                              paste0('S', formatC(abs(agb.lat), width = 2, flag = 0)), 
                              paste0('N', formatC(agb.lat, width = 2, flag = 0)))
        
        agb.fn <- NULL
        for (x in agb.lon.str) {
            for (y in agb.lat.str) agb.fn <- c(agb.fn, paste0(y, x, '_agb.tif'))
        }   # end for xy
        agb.files <- file.path(agb.path, agb.fn)

        # ------------------------------------------------------------------------ #
        ### read file as raster
        if (FALSE %in% file.exists(agb.files)) {
            missing.files <- basename(agb.files[!file.exists(agb.files)])
            stop(paste('grab.agb.GlobBiomass(): Missing AGB file in', agb.path, 
                       '\nPlease download files of "',  paste(missing.files, collapse = ', ')
                       , '" from http://globbiomass.org/wp-content/uploads/GB_Maps/Globbiomass_global_dataset.html...\n'))
        }   # end if
        
        cat(paste('grab.agb.GlobBiomass(): need to work on', length(agb.files), 
                  'file(s)...\n'))

        # deal with multiple files for GlobBiomass, DW, 04/19/2019
        if (length(agb.files) > 1) {
            
            cat('grab.agb.GlobBiomass(): loading and cropping multiple 100m AGB files...\n\n')
            agb.list <- list()
            for (b in 1 : length(agb.files)) 
                agb.list[[b]] <- crop(raster(agb.files[b]), ext)
            
            # merge multiple rasters
            cat('grab.agb.GlobBiomass(): merging multiple 100m AGB raster layers, please wait\n\n')
            agb.merge <- do.call(merge, agb.list); gc()

            # save it as tif 
            cat('grab.agb.GlobBiomass(): saving reformed AGB as tif file\n')
            writeRaster(agb.merge, out.fn, overwrite = T, format = 'GTiff', 
                        varname = 'AGB', xname = 'lon', yname = 'lat', 
                        longname = '100m AGB from GlobBiomass')

        } else {
            
            # crop spatial domain
            cat('grab.agb.GlobBiomass(): cropping 100m AGB...\n')
            sel.agb.rt <- crop(raster(agb.files[1]), ext)
            
            # save it as tif 
            writeRaster(sel.agb.rt, out.fn, overwrite = T, format = 'GTiff', 
                        varname = 'AGB', xname = 'lon', yname = 'lat', 
                        longname = '100m AGB from GlobBiomass')
                        
        }  # end if length(agb.files)
    }   # end if file.exists()

    return(out.fn)
}  # end of subroutine
