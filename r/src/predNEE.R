#' predNEE: script to downscale daily GPP and RECO to hrly NEE
#' this subroutine stores hourly NEE in a single day as tif files
#' @author: Dien Wu, 07/17/2019 

#' @param reg.name reg.name name, e.g., 'SaltLakeCity', without any space
#' @param reg.path path that stores GPP and NEE
#' @param yyyymm YYYYMM year + month

#' @param TA.path path for hourly air temperature fields
#' @param TA.field which Tair fields to use, only 'ERA5';
#'                   if using other fields, please write your own function
#' @param TA.varname the common characters of filenames in @param TA.path, 
#'                     e.g., TA.varname = '2T_' for '2T_201401.nc'

#' @param SSRD.path path for hourly surface shortwave radiation downward (SSRD)
#' @param SSRD.field which ssrd fields to use, only 'ERA5'
#' @param SSRD.varname the common characters of filenames in @param TA.path, 
#'                     e.g., SSRD.varname = 'SSRD_' for 'SSRD_201401.nc'
#' @param smurf_wd path for this bio model repository

# ---------------------------------------------------------------------------- #
predNEE <- function(reg.name = 'westernCONUS', 
                    reg.path, 
                    reco.dir = 'daily_mean_Reco_era5_era5',  
                    yyyymm = 201801,  # can be a vector, each core works on one
                    TA.path, 
                    TA.field = 'ERA5', 
                    TA.varname = '2T', 
                    SSRD.path, 
                    SSRD.field = c('ERA5', 'EPIC')[1], 
                    SSRD.varname = c('SSRD', NA)[1], 
                    smurf_wd = getwd()) {

  try({

    # Ensure dependencies are loaded for current node/process
    setwd(smurf_wd); source(file.path(smurf_wd, 'r/dependencies.r'), local = T)
    #if (!is.null(tmpdir)) raster::rasterOptions(tmpdir = tmpdir)

    # load all nc files that match 'yr' and look for GPP and Reco files
    yr <- substr(yyyymm, 1, 4)
    gpp.file <- list.files(reg.path, 'daily_mean_SIF_GPP_uncert', full.names = T)
    gpp.file <- gpp.file[grepl(yr, gpp.file)]

    # check whether files exist
    if (length(gpp.file) == 0) 
      stop('predNEE(): NO simulated SIF-based GPP found, see `main_script_GPP.r`\n')

    # stack all daily mean Reco files (in the same year) into one nc file
    # this process takes a while, will returns one nc file name
    # will be stored in reg.path
    cat('predNEE: Stacking all daily Reco into one nc file...\n')
    num <- as.numeric(lubridate::days_in_month(as.POSIXct(paste0(yyyymm, '01'), 
                                               format = '%Y%m%d', tz = 'UTC')))

    reco.file <- stack.reco(reg.name, reg.path, reco.dir, yyyymm, 
                            reco.file.num = num, overwriteTF = T)    
    
    reco.stk <- stack(reco.file, varname = 'Reco_mean')
    all.dates <- as.POSIXct(as.numeric(gsub('X', '', names(reco.stk))), 
                            origin = '1970-01-01 00:00:00', tz = 'UTC')
    all.timestr <- paste0(format(all.dates, format = '%Y%m%d'), '00')

    nee.path <- file.path(reg.path, paste0('hourly_flux_', tolower(SSRD.field)))
    dir.create(nee.path, showWarnings = F, recursive = T)

    # ------------------------------------------------------------------------ #
    ## get timestr in form of YYYYMMDDHH  
    cat(paste('\n\n# --- Working on YYYYMM:', yyyymm, 'for', reg.name, '--- #\n'))
    month <- formatC(substr(yyyymm, 5, 6), flag = 0, width = 2)
    mon.timestr <- all.timestr[substr(all.timestr, 5, 6) == month]

    # store hourly fluxes into nc files by months
    mean.gpp.stk <- mean.reco.stk <- mean.nee.stk <- NULL  
    for (tt in 1 : length(mon.timestr)) {

      cat('# ----------------------------------------------------- #')

      # call downscale.nee() for hourly downscaling
      # return list of hourly GPP, Reco and NEE CO2 fluxes
      hrly.list <- downscale.nee.hrly(mon.timestr[tt], gpp.file, reco.file, 
                                      TA.path, TA.field, TA.varname, SSRD.path, 
                                      SSRD.field, SSRD.varname)

      # store all hourly fluxes
      if (tt == 1) { 
          mean.gpp.stk  <- hrly.list$hrly_GPP_mean
          mean.reco.stk <- hrly.list$hrly_Reco_mean
          mean.nee.stk  <- hrly.list$hrly_NEE_mean

      } else {
          mean.gpp.stk  <- stack(mean.gpp.stk, hrly.list$hrly_GPP_mean)
          mean.reco.stk <- stack(mean.reco.stk, hrly.list$hrly_Reco_mean)
          mean.nee.stk  <- stack(mean.nee.stk, hrly.list$hrly_NEE_mean)
      }   # end if 

      gc()
    }   # end for tt


    # -------------------------- STORE hourly fluxes ----------------------- #
    cat(paste('\n\nStoring hourly fluxes as nc for', yyyymm, '\n'))
    nee.fn <- file.path(nee.path, paste0('hrly_mean_GPP_Reco_NEE_', reg.name, 
                                         '_', yyyymm, '.nc'))
    
    varnames  <- c('GPP_mean', 'Reco_mean', 'NEE_mean')
    varunits  <- rep('umol m-2 s-1', 3)
    longnames <- c('Hourly Mean Gross Primary Production (best estimates)', 
                   'Hourly Mean Ecosystem Respiration (best estimates)', 
                   'Hourly Mean Net Ecosystem Exchanges (best estimates)')

    # time format of names(mean.nee.stk), accuracy up to second in this case
    zformat  <- 'X%Y.%m.%d.%H.%M.%S' 
    stk.list <- list(mean.gpp.stk, mean.reco.stk, mean.nee.stk)
    names(stk.list) <- varnames

    # call save.raster2nc for storing multiple rasterStacks into one nc file
    save.raster2nc(varnames, varunits, longnames, zformat, stk.list, filename = nee.fn)
    
    return(nee.fn)  # return filename 
  })  # end of try()
}   # end of subroutine






# for sanity check, plot those scaling factors 
if (F) {

    selcol <- rev(c('#7F3B08', '#B35806', '#E08214', '#FDB863', '#FEE0B6', '#F7F7F7', 
                    '#D9F0D3', '#A6DBA0', '#5AAE61', '#1B7837', '#00441B'))    
    neeTheme <- rasterTheme(region = selcol)

    plot.nee <- hrly.list$hrly_NEE_mean; plot.nee[plot.nee == 0] <- NA
    n1 <- levelplot(plot.nee, layout = c(4, 6), at = seq(-55, 55, 5), 
                    par.settings = neeTheme, xlab = 'LONGITUDE', ylab = 'LATITUDE', 
                    max.pixels = 8e6, names.attr = as.character(seq(0, 23)))

    ggsave(arrangeGrob(n1), filename = '../paper3/gmd2020/test_era5_nee_20180701.png',
           width = 6, height = 8)

}


