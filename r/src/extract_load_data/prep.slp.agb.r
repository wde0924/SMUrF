#' subroutine to load or generate 500m GPP-SIF slopes and AGB data
#' before any urban gap-fills 
#' @author Dien Wu, 09/10/2019 

# separate C3 from C4, DW, 03/29/2020
prep.slp.agb <- function(prep.path, yr, slp.file, agb.file, ratio.file, 
                         minlon, maxlon, minlat, maxlat, lc.rt, overwriteTF = F) {

    reg.ext <- raster::extent(minlon, maxlon, minlat, maxlat)
    ext.str <- paste(minlon, maxlon, minlat, maxlat, sep = '_')

    # tif file names for gridded GPP-SIF slopes and 500 m AGB
    prep.fn <- file.path(prep.path, paste0('prep_gpp_', ext.str, '_', yr, '.nc'))

    # ------------------------------------------------------------------------ #
    ### need total three tif files, mean SLP, AGB
    # both at 500 m grid spacing that matches MCF
    if ( file.exists(prep.fn) & overwriteTF == F ) {
        cat(paste('prep.slp.agb(): Reading 500 m GPP-SIF slopes and AGB for region:', 
                   ext.str, '...\n'))

        slp.rt <- raster(prep.fn, varname = 'SLP')
        int.rt <- raster(prep.fn, varname = 'INT')
        cv.rt  <- raster(prep.fn, varname = 'GPP_CV')
        agb.rt <- raster(prep.fn, varname = 'AGB')

    } else {

        # load C4 fractions
        c4.rt <- raster(ratio.file); c3.rt <- 1 - c4.rt
        names(c3.rt) <- gsub('C4', 'C3', names(c4.rt))

        # read IGBP for a specific year
        igbp.stk <- stack(crop(lc.rt, reg.ext), c4.rt, c3.rt)
        igbp.df <- raster::as.data.frame(igbp.stk, xy = T)
        colnames(igbp.df) <- list('x', 'y', 'val', 'c4.frac', 'c3.frac')

        # get GPP-SIF slopes 
        raw.slp <- get.raw.slp(slp.file)
        slp.cro.c4 <- raw.slp[raw.slp$name == 'CRO_C4', ]
        slp.cro.c3 <- raw.slp[raw.slp$name == 'CRO_C3', ]

        # *** perform C3-C4 partitioning, val = 12 or 14 for CRO or CRO/NVM
        # DW, 03/29/2020
        cat('prep.slp.agb(): Calculating C3-C4 weighted mean GPP-SIF relation...\n')
        slp.all.df <- igbp.df %>% filter(!is.na(val)) %>%   # remove NA IGBP
                      left_join(raw.slp, by = 'val') %>% 
                      dplyr::select(-c('col', 'abbr'))
        
        # calculate the weighted mean slopes for CRO 
        rev.df <- slp.all.df %>% 
                  mutate(croTF = grepl('CRO', name), 
                         SLP = ifelse(croTF, c4.frac * slp.cro.c4$SLP + 
                                             c3.frac * slp.cro.c3$SLP, SLP), 

                         INT = ifelse(croTF, c4.frac * slp.cro.c4$INT + 
                                             c3.frac * slp.cro.c3$INT, INT), 

                         cv = ifelse(croTF, sqrt( (c4.frac * slp.cro.c4$cv) ^ 2 + 
                                                  (c3.frac * slp.cro.c3$cv) ^ 2 ), cv)
                        )   # end of mutate
        
        # ------------------------------------------------------------------------ #
        cat('prep.slp.agb(): Converting data frame back to rasterStack...\n')
        slp.stk <- rasterFromXYZ(rev.df %>% dplyr::select('x', 'y', 'SLP', 'INT', 'cv'))
        crs(slp.stk) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'

        # read biomass data and reproject 100m biomass to 500m
        cat('prep.slp.agb(): reprojecting 100m AGB to 500m Land cover grids...\n')
        agb.rt.100 <- crop(raster(agb.file), reg.ext) 
        agb.rt <- raster::projectRaster(agb.rt.100, slp.stk$SLP)    


        # -------------- save slopes, ints GPP cv, and AGB to nc file ----------
        varnames <- c('SLP', 'INT', 'GPP_CV', 'AGB')
        stk.list <- c(as.list(slp.stk), agb.rt)
        names(stk.list) <- varnames
        varunits <- c('(umol m-2 s-1) / (mW m-2 nm-1 sr-1)', 
                      '(umol m-2 s-1) / (mW m-2 nm-1 sr-1)', 
                      'unitless', 'tonne-C/ha')

        longnames <- c('GPP-CSIF slopes', 'GPP-CSIF intercepts', 
                       'GPP biome-specific fractional uncertainty (CV) based on FLUXNET GPP', 
                       'Reprojected 500m AGB from 100m GlobBiomass')

        prep.fn <- file.path(prep.path, paste0('prep_gpp_', ext.str, '_', yr, '.nc'))
        save.raster2nc(varnames, varunits, longnames, zformat = NA, stk.list, 
                       filename = prep.fn)
        
        slp.rt <- slp.stk$SLP
        int.rt <- slp.stk$INT 
        cv.rt  <- slp.stk$cv 
    }   # end if file.exists

    prep.stk <- stack(slp.rt, int.rt, cv.rt, agb.rt)
    names(prep.stk) <- c('SLP', 'INT', 'GPP_CV', 'AGB')

    cat('prep.slp.agb(): END of this subroutine...\n')
    return(prep.stk)
}   # end of subroutine