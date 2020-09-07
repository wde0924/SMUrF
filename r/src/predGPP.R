#' subroutine to estimate GPP
#' @author Dien Wu, 07/03/2019, latest modification on 03/28/2020

# ---------------------------------------------------------------------------- #
#' @param reg.name character string without any space for region name, e.g., 'westernCONUS'
#' @param minlon @param maxlon @param minlat @param maxlat spatial domain for Reco calculations
#' @param yr one numeric value (YYYY) indicating which year to work on 

#' @param lc.path input path for pre-processed MCD12Q1 files
#' @param lc.pattern unique character string of the MCD12Q1 file name
#' @param lc.max.yr latest year available for MCD12Q1 files, e.g., 2018
#' @param lc.res horizontal grid spacing of MCD12Q1 files
#' @param agb.file tif file for AGB, pre-processed by grab.agb.GlobBiomass()
#' @param ratio.file tif file for C4 ratio, pre-processed by prep.c4.ratio()

# ---------------------------------------------------------------------------- #
#' @param sif.prod character string indicating the name of adopted spatial 
#'                 SIF product, e.g., 'CSIFclear'
#' @param sif.path input path for this SIF product
#' @param sif.var variables to grab from SIF files, e.g., 'clear_daily_SIF'
#' @param sif.nd temporal resolution in days of this spatial SIF product
#'               e.g., 4 for 4 days
#' @param sif.res spatial grid spacing in degrees of this SIF product
#' @param sif.rmTF logical flag, TRUE for zero-ing negative SIF values

#' @param slp.file txtfile for GPP-SIF slopes, default is provided in the repository
#' @param slp.overwriteTF logical flag, whether to overwrite existing GPP-SIF 
#'                        slopes for natural biomes before urban gap-fill

# ---------------------------------------------------------------------------- #
#' @param gpp.path output path for storing modeled GPP fluxes
#' @param smurf_wd
#' @param tmpdir

# --------------------------------------------------------------------------- #
predGPP <- function(reg.name = 'westernCONUS',  # character
                    minlon = -125,       # in deg, negative for west hemisphere
                    maxlon = -95,  
                    minlat = 30, 
                    maxlat = 50, 
                    yr = 2018,  
                    
                    lc.path, 
                    lc.pattern = 'CONUS.tif', 
                    lc.max.yr = 2018, 
                    lc.res = 1/240,   # in deg
                    agb.file, 
                    ratio.file,     

                    sif.prod = 'CSIFclear', 
                    sif.path, 
                    sif.var = 'clear_daily_SIF',
                    sif.nd = 4,             # in days 
                    sif.res = 0.05,         # in deg
                    sif.rmTF = T,    # whether to force negative SIF as ZERO

                    gpp.path, 
                    smurf_wd,    # path for bio model repo
                    slp.file = file.path(smurf_wd, 'data/GPP-CSIF_stat_Wu.txt'), 
                    slp.overwriteTF = F, 
                    tmpdir = NA # temporary directory for processing large raster
                    ){


  try({

    # Ensure dependencies are loaded for current node/process
    setwd(smurf_wd); source(file.path(smurf_wd, 'r/dependencies.r'), local = T)
    if (!is.na(tmpdir)) raster::rasterOptions(tmpdir = tmpdir)
    reg.ext <- raster::extent(minlon, maxlon, minlat, maxlat)


    # -------------------------------------------------------------------------
    # 1) estimate GPP-SIF slopes with urban gap-fill 
    # -------------------------------------------------------------------------
    cat(paste('\npredGPP(): # --- STEP 1: Estimate GPP-SIF slopes for', yr, '--- #\n'))

    # prepare MCD12 IGBP land cover
    lc.rt <- prep.mcd12(lc.path, lc.pattern, yr, lc.max.yr, reg.name, reg.ext)

    # bias.corrTF logical flag, TRUE for performing urban bias correction
    #             according to Zhang et al. (2018), urban SIF may have an 
    #             underestimation of about 14.5%, thus we simply scale up 
    #             GPP-SIF slopes for MODIS-based urban areas
    bias.corrTF <- FALSE; if (grepl('CSIF', sif.prod)) bias.corrTF <- TRUE

    # Gap filling for urban SLOPEs by calling pred.slp.urban(), return a data.frame
    # see Wu et al. GMD, 2020 & calc.tree.type.frac.r for details, 03/28/2020
    slp.list <- pred.slp(gpp.path, reg.name, yr, minlon, maxlon, minlat, maxlat, 
                         slp.file, slp.overwriteTF, agb.file, ratio.file, 
                         lc.rt, lc.res, sif.res, bias.corrTF, smurf_wd) 
                # return a list of two rasterStacks

    prep.init   <- slp.list$stk.init    # initial slopes at 500 m
    prep.coarse <- slp.list$stk.coarse  # aggregate slopes at 0.05 deg
    cat('Done gap filling for GPP-SIF slopes of target city...\n\n')


    # -------------------------------------------------------------------------
    # 2) Generate 4-day mean GPP 
    # -------------------------------------------------------------------------
    cat(paste('predGPP(): # --- STEP 2: Generate GPP for every', sif.nd, 'day --- #\n'))

    ## Generate n day in a yr for temporal invest, format YYYYMMDDHH, in UTC
    all.timestr <- seq(as.Date(paste0(yr, '-01-01')), as.Date(paste0(yr, '-12-31')), sif.nd)
    all.timestr <- paste0(format(all.timestr, '%Y%m%d'), '00')
    all.date    <- as.POSIXct(all.timestr, format = '%Y%m%d', tz = 'UTC')

    ## loop over every 4 days in a particular yr
    gpp.mean.stk <- gpp.sd.stk <- sif.stk <- NULL    # initialize
    for (tt in 1 : length(all.timestr)) {

        timestr <- all.timestr[tt]
        if (tt %% 5 == 0) cat(paste('\n# ---- Working on date:', timestr, ';', 
                                    signif(tt / length(all.timestr)) * 100, 
                                    '% done --- #\n'))

        # ------------- 2.1 Grab spatial SIF and compute uncertainty --------- #
        ## grab two spatial SIF, they can be negative
        if (grepl('CSIF', sif.prod)) 
            sif.rt <- grab.csif(sif.path, timestr, ext = reg.ext, var = sif.var) 
        
        if (is.null(sif.rt)) stop(paste('predGPP(): No SIF file found for', 
                                  substr(timestr, 1, 8), 'Please check...\n'))
        if (sif.rmTF) sif.rt[sif.rt < 0] <- 0     # force negative CSIF to zero

        # ---------------------- 2.2 Compute gridded GPP -------------------- #
        # compute GPP with unit conversion to umol/m2/s, by calling compute.gpp()
        gpp.stk <- compute.gpp(sif.rt, prep.coarse)

        # save rasterLayer
        if (tt == 1) { 
            sif.stk      <- sif.rt 
            gpp.mean.stk <- gpp.stk$GPP_mean
            gpp.sd.stk   <- gpp.stk$GPP_sd 

        } else { 
            sif.stk      <- stack(sif.stk,      sif.rt)
            gpp.mean.stk <- stack(gpp.mean.stk, gpp.stk$GPP_mean) 
            gpp.sd.stk   <- stack(gpp.sd.stk,   gpp.stk$GPP_sd) 
        }   # end if tt == 1

        # generate a plot for one summertime day
        if (tt == 50) {     
            plot.stk <- gpp.stk
            g1 <- gridExtra::grid.arrange(levelplot(plot.stk, maxpixel = 1e6))
            title <- paste('4-day mean 0.05◦ GPP and uncertainty [umol m-2 s-1] for', 
                           reg.name, 'on', timestr)

            g1 <- annotate_figure(g1, top = title)
            fn <- file.path(gpp.path, paste0('GPP_SD_', reg.name, '_', timestr, '.png'))
            ggsave(g1, filename = fn, width = 10, height = 6)
        }   # end if tt == 50

    }   # end for tt for all days

    # -------------------------------------------------------------------------
    # 3) save GPP in nc files 
    # -------------------------------------------------------------------------
    cat('predGPP(): # --- STEP 3: Saving GPP in nc files --- #\n')
    gpp.fn <- file.path(gpp.path, paste0('daily_mean_SIF_GPP_uncert_', reg.name, 
                                         '_', yr, '.nc'))

    varnames  <- c('SIF_mean', 'GPP_mean', 'GPP_sd')
    varunits  <- c('mW m−2 nm−1 sr−1', rep('umol m-2 s-1', 2) )
    longnames <- c('Daily Mean clear-sky CSIF [Zhang et al., 2018]', 
                   'Daily Mean Gross Primary Productions (best estimates)', 
                   '1-sigma Uncertainty of Daily Mean Gross Primary Productions')
    zformat <- 'X%Y.%m.%d'

    # assign correct layer names 
    names(sif.stk) <- names(gpp.mean.stk) <- names(gpp.sd.stk) <- all.date
    
    # order of this list should match all variables above e.g., varnames
    stk.list <- list(sif.stk, gpp.mean.stk, gpp.sd.stk)
    names(stk.list) <- varnames

    # call save.raster2nc for storing multiple rasterStacks into one nc file
    save.raster2nc(varnames, varunits, longnames, zformat, stk.list, filename = gpp.fn)

  })  # end of try()

}   

# end of subroutine




# ---------------------------------------------------------------------------- #
# sanity check 
if (F) {
    
    plot.gpp <- gpp.stk$GPP_mean
    library(gridExtra)
    col <- rasterTheme(region = rev(terrain.colors(12)))
    #e1 <- extent(-112.3, -111.4, 40, 41); e2 <- extent(-119.5, -117, 32, 35); e3 <- extent(-123, - 121.5, 46.8, 48.5)
    #lab <- c('SLC', 'LA', 'SEA')
    e1 <- extent(-87, -85, 39, 41); e2 <- extent(-76, -74, 39, 41); e3 <- extent(-72, -70, 41, 43)
    lab <- c('Indy', 'Baltimore-DC', 'Boston')

    s1 <- grid.arrange(levelplot(crop(plot.gpp, e1), at = seq(-1, 30, 1), par.settings = col))
    s2 <- grid.arrange(levelplot(crop(plot.gpp, e2), at = seq(-1, 30, 1), par.settings = col))
    s3 <- grid.arrange(levelplot(crop(plot.gpp, e3), at = seq(-1, 30, 1), par.settings = col))

    ss <- ggarrange(s1, s2, s3, ncol = 3, widths = c(1.03, 1.03, 1.01), labels = lab) 
    ss <- annotate_figure(ss, top = 'GPP with urban gap-fills zoomed into selected cities')
    ggsave(ss, filename = '../GPP_zoom_east.png', width = 13, height = 6)

}   # end if

