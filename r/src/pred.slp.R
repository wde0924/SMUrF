#' script to gap fill slopes for urban
#' @author Dien Wu, 04/19/2019, last modification on 03/28/2020

#' see @function predGPP.R for explanation of input parameters

# half OSHR and half GRASS
pred.slp <- function(gpp.path, 
                     reg.name, 
                     yr, 
                     minlon, 
                     maxlon, 
                     minlat, 
                     maxlat, 

                     slp.file, 
                     slp.overwriteTF = F, 
                     agb.file, 
                     ratio.file, 

                     lc.rt, 
                     lc.res = 1/240, 
                     sif.res = 0.05, 
                     bias.corrTF = T, 
                     smurf_wd) {

    # ----------------------------------------------------------------------- #
    # call prep.slp.agb to load or generate 500m GPP-SIF slopes and AGB based on MODIS
    # will return a list of one rasterStack and one rasterLayer
    # separate C3 from C4, DW, 03/29/2020
    prep.stk <- prep.slp.agb(gpp.path, yr, slp.file, agb.file, ratio.file, 
                             minlon, maxlon, minlat, maxlat, lc.rt, slp.overwriteTF) 
    
    # get empirical relation and add a row for nearfield grids if IGBP = 0
    raw.slp <- get.raw.slp(slp.file)
   
    # add 500 m land cover IGBP values and convert rasterStack to data frame
    input.stk <- stack(lc.rt, prep.stk)
    names(input.stk) <- list('val', 'SLP', 'INT', 'GPP_CV', 'AGB')
    input.df <- raster::as.data.frame(input.stk, xy = T) %>% 
                left_join(raw.slp[, c('val', 'abbr')], by = 'val')
    

    # ------------------------------------------------------------------------ #
    # we treat cities as combinations of trees and non-trees (e.g., 50% OSHR and 50% GRASS)
    # 1) calculate tree fractions based on AGB
    # ------------------------------------------------------------------------ #
    cat('pred.slp(): 1) estimate and assign tree fractions to MODIS-based urban areas...\n')
    urban.df <- input.df %>% filter(abbr == 'URB') %>% 
                mutate(tot.tree.frac = (1.68 * AGB ^ 0.577 + 33.5) / 100, 
                       tot.tree.frac = ifelse(tot.tree.frac > 1, 1, tot.tree.frac),  
                       OSHR.frac = (1 - tot.tree.frac) / 2, 
                       GRA.frac = (1 - tot.tree.frac) / 2)  # fraction in %


    # ------------------------------------------------------------------------ #
    # 2) call calc.tree.type.frac() to figure out tree types including ENF, DBF, & EBF
    # ------------------------------------------------------------------------ #
    # first select urban grids
    cat('pred.slp(): 2) Assigning WEIGHTED GPP-SIF slopes to MODIS-derived URBAN...\n')
    urban.tree.type <- get.tree.type.frac(urb.lat = urban.df$y, smurf_wd) 
    urban.frac.df <- cbind(urban.df, urban.tree.type) %>%
                     mutate(DBF.frac = DBF_rel * tot.tree.frac, 
                            DNF.frac = DNF_rel * tot.tree.frac, 
                            EBF.frac = EBF_rel * tot.tree.frac,
                            ENF.frac = ENF_rel * tot.tree.frac, 
                            MF.frac  = MF_rel * tot.tree.frac) %>% 
                     dplyr::select(-c(contains('_rel'), 'match.lat', 'urb.lat')) 
        
    # reframe the data frame of slopes and intercepts 
    reframe.slp <- function(name = 'DBF', raw.slp) {
        rev.slp <- raw.slp[raw.slp$name == name, ]
        rev.slp <- data.frame(slp = rev.slp$SLP, int = rev.slp$INT, 
                              cv = rev.slp$cv, igbp = rev.slp$val)
        colnames(rev.slp) <- paste0(name, '.', colnames(rev.slp))
        return(rev.slp)
    }   # end of function 

    mf.slp <- reframe.slp('MF', raw.slp)
    dbf.slp <- reframe.slp('DBF', raw.slp); dnf.slp <- reframe.slp('DNF', raw.slp)
    ebf.slp <- reframe.slp('EBF', raw.slp); enf.slp <- reframe.slp('ENF', raw.slp)
    oshr.slp <- reframe.slp('OSHR', raw.slp); gras.slp <- reframe.slp('GRA', raw.slp)

    # calculate the biome-weighted SLP and INT
    urban.slp.df <- cbind(urban.frac.df, mf.slp, dbf.slp, dnf.slp, ebf.slp, 
                                         enf.slp, oshr.slp, gras.slp) %>% 

                    mutate(SLP.wgt = DNF.frac * DNF.slp + DBF.frac * DBF.slp + 
                                     EBF.frac * EBF.slp + ENF.frac * ENF.slp + 
                                     MF.frac * MF.slp + OSHR.frac * OSHR.slp + 
                                     GRA.frac * GRA.slp, 

                           INT.wgt = DNF.frac * DNF.int + DBF.frac * DBF.int + 
                                     EBF.frac * EBF.int + ENF.frac * ENF.int + 
                                     MF.frac * MF.int + OSHR.frac * OSHR.int +
                                     GRA.frac * GRA.int, 

                           GPP_CV.wgt = sqrt( (DNF.frac * DNF.cv) ^ 2 + 
                                              (DBF.frac * DBF.cv) ^ 2 + 
                                              (EBF.frac * EBF.cv) ^ 2 + 
                                              (ENF.frac * ENF.cv) ^ 2 + 
                                              (MF.frac  * MF.cv)  ^ 2 + 
                                              (OSHR.frac * OSHR.cv) ^ 2 + 
                                              (GRA.frac  * GRA.cv) ^ 2
                                            ) ) %>% 

                    dplyr::select(x, y, ends_with('.wgt'))
    colnames(urban.slp.df) <- gsub('.wgt', '', colnames(urban.slp.df))

    # bias correction for urban grids, since CSIF may be underestimated by 14.5%
    if (bias.corrTF) urban.slp.df <- urban.slp.df %>% mutate(SLP = SLP * 1.145)          


    # ------------------------------------------------------------------------ #
    # 3) combine slopes
    # ------------------------------------------------------------------------ #
    # convert back to raster and merge with slopes for natural biomes
    # and fix the extend by calling raster::extend() 
    cat('pred.slp(): 3) Merging Urban slopes and fractions for cities back to initial ones...\n')
    urban.slp.brk <- rasterFromXYZ(urban.slp.df[, c('x', 'y', 'SLP', 'INT', 'GPP_CV')])
    crs(urban.slp.brk) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    urban.slp.brk <- raster::extend(urban.slp.brk, extent(prep.stk))
    
    # merge gap-filled rasters back to initial rasters
    merge.slp <- raster::merge(urban.slp.brk$SLP, prep.stk$SLP)
    merge.int <- raster::merge(urban.slp.brk$INT, prep.stk$INT)
    merge.cv  <- raster::merge(urban.slp.brk$GPP_CV, prep.stk$GPP_CV)
    stk.init <- stack(merge.slp, merge.int, merge.cv)
    names(stk.init) <- c('SLP', 'INT', 'GPP_CV')


    # ------------------------------------------------------------------------ #
    # 4) *** save 500m land cover fractions (OSHR, GRAS, ENF, EBF, DBF, DNF, MF)
    # and store GPP-CSIF slopes before and after urban gap-fills 
    # ------------------------------------------------------------------------ #
    cat('pred.slp(): 4) Saving urban veg fractions and GPP-SIF relation in nc file...\n')
    
    # convert gap-filled veg fractions from data fraom to rasterStacks
    # and fix rasterBrick's extent
    sel.frac.df <- urban.frac.df %>% dplyr::select(c(x, y, contains('.frac')))
    urban.frac.brk <- rasterFromXYZ(sel.frac.df)
    crs(urban.frac.brk) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    names(urban.frac.brk) <- gsub('.frac', '_frac', colnames(sel.frac.df)[-c(1, 2)])
    urban.frac.brk <- raster::extend(urban.frac.brk[[-1]], extent(prep.stk))

    # ------ 
    output.stk <- stack(merge.slp, merge.int, urban.frac.brk)
    crs(output.stk) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    names(output.stk)[1:2] <- c('SLP_gapfill', 'INT_gapfill')

    output.fn <- file.path(gpp.path, paste0('veg_fraction_slp_', reg.name, '_', yr, '.nc'))
    varnames <- names(output.stk)
    stk.list <- as.list(output.stk)
    names(stk.list) <- varnames
    varunits <- rep('unitless', nlayers(output.stk))
    longnames <- c('GPP-CSIF slopes after urban gap-fill', 
                   'GPP-CSIF intercept after urban gap-fill', 
                   'Urban OSHR fraction', 'Urban GRA fraction', 
                   'Urban DBF fraction', 'Urban DNF fraction', 
                   'Urban EBF fraction', 'Urban ENF fraction', 
                   'Urban MF fraction')

    # call save.raster2nc for storing multiple rasterStacks into one nc file
    #' @param stk.list needs to have a list name
    save.raster2nc(varnames, varunits, longnames, zformat = NA, stk.list, 
                   filename = output.fn)


    # ------------------------------------------------------------------------ #
    # 5) aggregate slopes + uncertainty (calculate RMSE) to 0.05 deg, DW, 06/10/2019
    # ------------------------------------------------------------------------ #
    cat('pred.slp(): 5) Aggregating slopes from 500 m to SIF"s resolution...\n')
    coarse.slp <- raster::aggregate(merge.slp, fact = sif.res / lc.res)
    coarse.int <- raster::aggregate(merge.int, fact = sif.res / lc.res)
    coarse.cv <- sqrt(raster::aggregate(merge.cv ^ 2, fact = sif.res / lc.res))
    stk.coarse <- stack(coarse.slp, coarse.int, coarse.cv)
    names(stk.coarse) <- c('SLP', 'INT', 'GPP_CV')

    slp.list <- list(stk.init = stk.init, stk.coarse = stk.coarse)
    return(slp.list)

}   

# end of subroutine

