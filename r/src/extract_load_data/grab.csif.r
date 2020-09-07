# subroutine to grab CSIF from Zhang et al., 2018, BG
# by Dien Wu, 12/21/2018

# timestr in form of YYYYMMDDHH
# ext: generated from raster::extent(), e.g., minlon, maxlon, minlat, maxlat
grab.csif <- function(csif.path, timestr, ext = NULL, nhrs = NULL,
                      var = c('clear_inst_SIF', 'clear_daily_SIF', 
                              'all_inst_SIF', 'all_daily_SIF')[2], 
                      form = c('df', 'raster')[2]) {
    
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

    if (length(csif.file) == 0) {
        cat(paste('grab.csif(): No CSIF file found for', search.timestr, 
                  ', please check...\n')); return()

    } else {
      
      # grab and crop CSIF according to 'var'
      sel.csif <- crop(stack(csif.file, varname = var), ext)
      if (nlayers(sel.csif) > 1) sel.csif <- mean(sel.csif)

      if (form == 'df') {
        sel.csif <- raster::as.data.frame(sel.csif, xy = T)
        colnames(sel.csif) <- list('lon', 'lat', 'CSIF')
      } # end if form

      return(sel.csif)
    } # end if 

} # end of subroutine
