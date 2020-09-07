#' script to loop over and read each daily mean Reco, stack them together, and
#' store into one nc file, before downscaling GPP, Reco, NEE fluxes
#' @author Dien Wu, 08/13/2019

#' @param reg.name region name, e.g., 'westernCONUS'
#' @param yyyymm which year and month to work on 
#' @param reco.path path that stores daily reco files, see 'meain_script_Reco.r'
#' @param reco.file.num number of daily files that one should have
#' @param overwriteTF: whether to re-loading and storing all daily files 
#'                     to one nc file, even if already reco.file exists

stack.reco <- function(reg.name = 'westernCONUS', reg.path, reco.dir, 
                       yyyymm = 201807, reco.file.num = 31, overwriteTF = F) {
    
    # this is the path that is defined when running main_script_Reco.r
    reco.path <- file.path(reg.path, reco.dir, substr(yyyymm, 1, 4))
    fn <- file.path(reg.path, paste0('daily_mean_Reco_uncert_', reg.name, '_', yyyymm, '.nc'))
    
    if (!file.exists(fn) | overwriteTF == T) {

        cat('stack.reco(): NO stacked nc file found for Reco OR NEED to overwrite
             Start aggregating daily mean Reco in a year to one single nc file...\n')

        # look for all reco files
        reco.files <- list.files(reco.path, pattern = 'daily_mean_Reco_uncert', full.names = T)
        reco.files <- reco.files[grepl(yyyymm, basename(reco.files))]
        reco.files <- reco.files[basename(reco.files) != basename(fn)]

        if (length(reco.files) == reco.file.num) {  # if those files exist

            mean.stk <- raster(reco.files[1], varname = 'Reco_mean')
            sd.stk   <- raster(reco.files[1], varname = 'Reco_sd')
            z.values <- getZ(mean.stk)      # get Z values for timestr

            for (r in 2 : length(reco.files)) {
                if (r %% 10 == 0) cat(paste(signif(r / reco.file.num) * 100, '% done...\n'))
                mean.tmp <- raster(reco.files[r], varname = 'Reco_mean')
                mean.stk <- stack(mean.stk, mean.tmp)
                sd.stk   <- stack(sd.stk, raster(reco.files[r], varname = 'Reco_sd'))
                z.values <- c(z.values, getZ(mean.tmp))     
            }   # end for r

            z.times <- as.POSIXct(z.values, origin = '1970-01-01 00:00:00', tz = 'UTC')
            names(mean.stk) <- names(sd.stk) <- z.times

            # store rasterStacks into one nc file
            varnames  <- c('Reco_mean', 'Reco_sd')
            varunits  <- rep('umol m-2 s-1', 2)
            longnames <- c('Daily Mean Ecosystem Respiration (best estimates)', 
                           '1-sigma Uncertainty of Daily Mean Ecosystem Respiration')
            zformat <- 'X%Y.%m.%d'
            reco.list <- list(mean.stk, sd.stk)
            names(reco.list) <- varnames

            # save all daily reco into one nc file
            save.raster2nc(varnames, varunits, longnames, zformat, 
                           stk.list = reco.list, filename = fn)

        } else {
            stop(paste('stack.reco(): Cannot combine daily mean .nc files into one file yet\n',
                       '*** Need', reco.file.num, 'files for', yyyymm, '; you have', length(reco.files), 'files ***\n'))

        }   # end if length() 

    }   # end if stack files 

    return(fn)

}   # end of function 