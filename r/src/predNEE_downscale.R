#' predNEE: script to downscale daily GPP and RECO to hrly NEE
#' this subroutine stores hourly NEE in a single day as tif files
#' @author: Dien Wu, 07/17/2019 
#' @param reg.path path that stores GPP and NEE
#' 
#' @param smurf_wd path for this bio model repository
#' @param tmpdir temporary space for large raster calculations 

#' updates: 


# ---------------------------------------------------------------------------- #
predNEE.downscale <- function(site = 'SaltLakeCity',
                              reg.path,   
                              timestr = 20180716,  # only one YYYYMMDD

                              minlon = -115, 
                              maxlon = -109, 
                              minlat = 38, 
                              maxlat = 44,
                              
                              vcf.path, 
                              vcf.pattern, 
                              smurf_wd = getwd(),  # working directory 
                              tmpdir = NA) {

  try({

    # Ensure dependencies are loaded for current node/process
    setwd(smurf_wd); source(file.path(smurf_wd, 'r/dependencies.r'), local = T)
    if (!is.na(tmpdir)) raster::rasterOptions(tmpdir = tmpdir)

    flux.files <- list.files(reg.path, 'hrly_mean_GPP_Reco_NEE', full.names = T)
    if (length(flux.files) == 0) stop('predNEE.downscale(): no NEE files found for this region')

    flux.file <- flux.files[grepl(substr(timestr, 1, 6), flux.files)]
    if (length(flux.file) == 0) stop('predNEE.downscale(): no NEE files found for this region & timestr')

    #' call downscale.nee.spatial() to perform spatial downscaling using VCF 
    site.ext <- extent(minlon, maxlon, minlat, maxlat)
    site.path <- file.path(reg.path, site)
    dir.create(site.path, showWarnings = F)
    
    source('r/dependencies.r')  
    hrly.list <- downscale.nee.spatial(site, timestr, site.ext, flux.file, 
                                       vcf.path, vcf.pattern)
    print(length(hrly.list))

    if (F) {
      selcol <- rev(c('#7F3B08', '#B35806', '#E08214', '#FDB863', '#FEE0B6', '#F7F7F7', 
                      '#D9F0D3', '#A6DBA0', '#5AAE61', '#1B7837', '#00441B'))    
      colTheme <- rasterTheme(region = selcol)
      daily.ds <- stack(lapply(hrly.list, mean))
      d1 <- levelplot(crop(daily.ds[[3]], site.ext - 4), at = seq(-2.5, 2.5, length = 20), 
                      par.settings = colTheme, xlab = 'LONGITUDE', ylab = 'LATITUDE')

    }

    # -------------------------- STORE hourly fluxes ----------------------- #
    cat(paste('\n\npredNEE.downscale(): Storing hourly fluxes as nc for', timestr, '\n'))
    ds.path <- file.path(site.path, 'SMUrF_downscale'); dir.create(ds.path, showWarnings = F)
    nee.fn <- file.path(ds.path, paste0('hrly_fluxes_', site, '_', timestr, '.nc'))
    
    varunits  <- rep('umol m-2 s-1', 3)
    longnames <- c('VCF-downscaled hourly mean Gross Primary Production (best estimates)', 
                   'VCF-downscaled hourly mean Ecosystem Respiration (best estimates)', 
                   'VCF-downscaled hourly mean Net Ecosystem Exchanges (best estimates)')
    varnames <- names(hrly.list)

    # time format of names(mean.nee.stk), accuracy up to second in this case
    zformat  <- 'X%Y%m%d%H' 

    # call save.raster2nc for storing multiple rasterStacks into one nc file
    save.raster2nc(varnames, varunits, longnames, zformat, stk.list = hrly.list,
                   filename = nee.fn)

    return(nee.fn)  # return filename 
  })  # end of try()
}   # end of subroutine


