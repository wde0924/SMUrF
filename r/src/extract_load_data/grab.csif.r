# subroutine to grab CSIF from Zhang et al., 2018, BG
# by Dien Wu, 12/21/2018
#' updated by @author: Sabrina Madsen-Colford, 07/11/2022

# timestr in form of YYYYMMDDHH
# ext: generated from raster::extent(), e.g., minlon, maxlon, minlat, maxlat
grab.csif <- function(csif.path, timestr, sif.temp, TA.path, TA.varname,
                      ext = NULL, nhrs = NULL,
                      var = c('clear_inst_SIF', 'clear_daily_SIF', 
                              'all_inst_SIF', 'all_daily_SIF')[2], 
                      form = c('df', 'raster')[2], yr2=2018) {
    
    # call a funcion to get DOY from YYYYMMDD, given day intervals
    if (nchar(timestr) != 10) stop('grab.csif(): timestr needs to be in form of YYYYMMDDHH')
    time.info <- get.doy.from.timestr(timestr, nday = 4, nhrs = nhrs)
    search.timestr <- time.info$find.yrdoy

    
    # select the ones that need for scaling GPP 
    csif.files <- list.files(path = csif.path, pattern = '.nc', 
                             recursive = T, full.names = T)
    if (is.null(nhrs)) {
      csif.file <- csif.files[grep(search.timestr, csif.files)]

    } else {
      ref.csif.files <- gsub('.v2', '', csif.files)
      file.timestr <- strsplit.to.df(basename(ref.csif.files), sep = '\\.')$V5
      csif.file <- csif.files[file.timestr >= search.timestr & 
                              file.timestr <= time.info$nhrs.yrdoy]
    } # end if
    
    # SM: if there is a shoreline-corrected SIF use that, 10/23/2024
    if(length(csif.file)>1){
      csif.file <- csif.file[grepl('shore_weighted_corr',csif.file)]
      print(paste0('using shore corrected SIF: ',csif.file))
    }else{
      print(paste0('No shore correction available, ',csif.file))
    }
    
    if (length(csif.file) == 0) {
        cat(paste('grab.csif(): No SIF file found for', search.timestr, 
                  ', please check...\n')); return()

    } else {
      
      # grab and crop SIF according to 'var'
      sel.csif <- crop(stack(csif.file, varname = var), ext)
      if (nlayers(sel.csif) > 1) sel.csif <- mean(sel.csif)
      
      # SM: Also grab and crop temperature data 07/11/2022
      # Load hourly TA (deg C) and incoming SW from ERA5
      # use nhrs = 23 to grab all hourly variables in a day 
      
      # create 4day mean interval that matches GPP 
      yr <- substr(search.timestr, 1, 4)
      date4 <- seq(as.Date(paste0(yr, '-01-01')), as.Date(paste0(yr, '-12-31')), by = '4 days')
      timestr4 <- paste0(format(date4, format = '%Y%m%d'), '00')
      # find the start time of the 4-day interval that "timestr" fall into 
      find.timestr <- timestr4[findInterval(time.info$find.timestr, timestr4)]
      
      # SM: If it is outside of the growing season set SIF to 0, 10/07/2022
      # NOTE THIS WILL NEED TO BE CHANGED OUTSIDE OF THE GTA!!!
      if ((search.timestr-yr2*1000)<90 | (search.timestr-yr2*1000)>340){    
        sel.csif<-sel.csif*0
      }else if (sif.temp==TRUE) {     # hourly air temp in UTC
        TA.brk <- prep.era5(TA.path, TA.varname, timestr = find.timestr, 
                            ext, nhrs = 24 * 4)
        TA.pj   <- raster::projectRaster(TA.brk, sel.csif) # in degC
        
        TA.mean <- mean(TA.pj) #If mean air temperature is <0
        sel.csif[0 >= TA.mean]<-TA.mean[0 >= TA.mean]*0 #set SIF to 0
      }

      if (form == 'df') {
        sel.csif <- raster::as.data.frame(sel.csif, xy = T)
        colnames(sel.csif) <- list('lon', 'lat', 'CSIF')
      } # end if form

      return(sel.csif)
    } # end if 

} # end of subroutine
