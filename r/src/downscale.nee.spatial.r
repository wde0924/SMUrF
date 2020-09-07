#' subroutine to spatially downscale NEE using VCF from MOD44B
#' @author Dien Wu, 11/30/2019 

downscale.nee.spatial <- function(site, timestr, site.ext, flux.file, vcf.path, 
                                  vcf.pattern) {

    #' grab and crop hourly fluxes based on @param dlat and @param dlon 
    #' return 24 rasterLayers with raster names indicating YYYYMMDDHH
    cat('downscale.nee.spatial(): Reading 0.05deg bio fluxes...\n')
    gpp.stk <- crop.smurf.hrly.flux(timestr = paste0(timestr, '00'), flux.file, 
                                    site.ext, varname = 'GPP_mean', nhrs = 23)
    reco.stk <- crop.smurf.hrly.flux(timestr = paste0(timestr, '00'), flux.file, 
                                    site.ext, varname = 'Reco_mean', nhrs = 23)
    nee.stk <- reco.stk - gpp.stk

    # ------------------------------------------------------------------------ #
    # spatially downscaling using VCF
    # prepare gridded scaling factors, same for each same year
    # get 250m vegetated fractions
    yr <- substr(timestr, 1, 4)
    fveg.stk <- prep.mod44b(vcf.path, vcf.pattern, yr, reg.ext = site.ext)   # 250m 
    fveg.rt <- fveg.stk$veg_frac    # at 250m 

    # then reproject those fractions to get mean fractions at 1km
    fveg.pj <- projectRaster(fveg.rt, reco.stk[[1]], method = 'bilinear')  
    fveg.pj <- aggregate(fveg.rt, fact = 4)        # aggregate from 250m to 1km

    if (F) v1 <- levelplot(crop(fveg.pj, site.ext - 3), at = seq(0, 100, 5), 
                           main = paste('1km MODIS VCF for', site), 
                           par.settings = rasterTheme(region = rev(terrain.colors(20))), 
                           xlab = 'LONGITUDE', ylab = 'LATITUDE', margin = F)

    # ------------------------------------------------------------------------ #
    cat('downscale.nee.spatial(): Bilinear calculations to GPP, Reco, and fveg...it takes a while\n')
    # bilinear interpolate GPP, Reco, and veg fractions to 1km 
    gpp.bilinear  <- projectRaster(gpp.stk, fveg.pj)   
    reco.bilinear <- projectRaster(reco.stk, fveg.pj)

    fveg.bilinear <- projectRaster(aggregate(fveg.rt, fact = 24), fveg.pj) 
    fveg.bilinear[fveg.bilinear < 0] <- 0       # blurring one
    fveg.sf <- fveg.pj / fveg.bilinear; fveg.sf[fveg.sf == Inf] <- 0 

    # scale fluxes to 1km 
    cat('downscale.nee.spatial(): Downscaling bio fluxes to 1km...it takes a while\n')
    gpp.1km  <- gpp.bilinear * fveg.sf
    reco.1km <- reco.bilinear * fveg.sf

    # in order to store data in .nc, no allowance for infinite values
    gpp.1km[gpp.1km < 0] <- 0; gpp.1km[gpp.1km == -Inf] <- 0
    reco.1km[reco.1km < 0] <- 0; reco.1km[reco.1km == -Inf] <- 0 
    nee.1km <- reco.1km - gpp.1km
    names(gpp.1km) <- names(reco.1km) <- names(nee.1km) <- names(gpp.stk)

    # return all hourly GPP, Reco and NEE
    hrly.list <- list(gpp.1km, reco.1km, nee.1km)
    names(hrly.list) <- c('GPP_mean_1km', 'Reco_mean_1km', 'NEE_mean_1km')
    
    return(hrly.list)
}   

# end of function








######
if (F) {

    m1 <- ggplot.map(map = 'black', center.lat = (minlat + maxlat) / 2, 
                     center.lon = (minlon + maxlon) / 2, zoom = 9)
    g1 <- ggmap.gpp(gpp.stk[[20]], site, timestr, m1)
    g2 <- ggmap.gpp(gpp.bilinear[[20]], site, timestr, m1, res = '1km mean')
    g3 <- ggmap.gpp(gpp.1km[[20]], site, timestr, m1, res = '1km scaled')
    s1 <- ggmap.gpp(fveg.pj, site, timestr, m1, res = '1km fveg')
    gg <- ggarrange(g1, g2, s1, g3, ncol = 2, nrow = 2)
    ggsave(gg, filename = paste0('../downscale_GPP_fveg_', site, '_zoom.png'), 
                width = 10, height = 10)


    #### for plotting, just plot NEE at 20 UTC, on day 15 of each month
    hr <- 20; hr <- formatC(hr, width = 2, flag = 0)
    day <- 15; day <- formatC(day, width = 2, flag = 0)
    l <- which(names(gpp.stk) == paste0('X', timestr, hr))
    
    ###
    init <- stack(-gpp.stk[[l]], reco.stk[[l]], nee.stk[[l]])
    fine <- stack(-gpp.1km[[l]], reco.1km[[l]], nee.1km[[l]])
    sel.init <- crop(init, site.ext - 4); sel.init[sel.init == 0] <- NA
    sel.fine <- crop(fine, site.ext - 4); sel.fine[sel.fine == 0] <- NA
    names(sel.init) <- names(sel.fine) <- c('GEE', 'Reco', 'NEE')
    
    selcol <- rev(c('#7F3B08', '#B35806', '#E08214', '#FDB863', '#FEE0B6', '#F7F7F7', 
                '#D9F0D3', '#A6DBA0', '#5AAE61', '#1B7837', '#00441B'))    
    colTheme <- rasterTheme(region = selcol)

    zmax1 <- ceiling(max(abs(getValues(sel.init)), na.rm = T))
    zmax2 <- ceiling(max(abs(getValues(sel.fine)), na.rm = T))
    zmax1 <- zmax2 <- 14
    a1 <- levelplot(sel.init[[3]], at = seq(-zmax1, zmax1, length = 25), #layout = c(3, 1), 
                    par.settings = colTheme, xlab = 'LONGITUDE', ylab = 'LATITUDE', margin = F, 
                    main = paste('5km hourly fluxes for', site, 'on', timestr, hr, 'UTC'))
    a2 <- levelplot(sel.fine[[3]], at = seq(-zmax2, zmax2, length = 25), #layout = c(3, 1), 
                    par.settings = colTheme, xlab = 'LONGITUDE', ylab = 'LATITUDE', margin = F, 
                    main = paste('1km hourly fluxes using MODIS VCF for', site, 'on', timestr, hr, 'UTC'))
    aa <- arrangeGrob(a1, a2)  #;aav <- grid.arrange(aa, v1, ncol = 2)
    
    fn <- paste0('downscale_5vs1km_', site, '_', paste0(timestr, hr), '.png')
    ggsave(aa, filename = file.path(site.path, fn), width = 11, height = 10)

}