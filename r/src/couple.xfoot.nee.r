# subroutine to couple hourly xfootprint with NEE fluxes, DW, 11/21/2019

couple.xfoot.nee <- function(site, timestr, store.path, nee.path, foot.res, 
                             foot.pattern = '0.05x0.05_foot.nc', foot.nhrs = -24, 
                             nee.varname = 'NEE_mean', overwriteTF = F) {
    
    fn <- file.path(store.path, 'grd', paste0('foot_nee_', site, '_', timestr, 
                                       '_', foot.res, '.grd'))
    
    # load footprint files
    out.dirs   <- list.files(store.path, 'out_20', full.names = T)
    foot.path  <- out.dirs[grepl(timestr, out.dirs)]
    foot.files <- list.files(foot.path, foot.pattern, full.names = T, recursive = T)
    cat(paste('couple.xfoot.nee(): total', length(foot.files), 'footprint files\n'))

    foot.info <- strsplit.to.df(basename(foot.files))[, c(2, 3)] %>% rename(lon = V2, lat = V3)
    foot.ext <- extent(raster(foot.files[1]))

    # need to organize footprints by the order of their receptor latitudes
    order.indx <- order(foot.info$lat)
    foot.info  <- foot.info[order.indx,] %>% mutate_all(funs(as.numeric), colnames(foot.info)) 
    foot.files <- foot.files[order.indx]
    recp.lat   <- as.numeric(foot.info$lat)
    recp.lon   <- as.numeric(foot.info$lon)

    if (file.exists(fn) && overwriteTF == FALSE) {
        cat('couple.xfoot.nee(): found file, calculate the total XCO2.bio\n')
        sum.fn.stk <- stack(fn)
        names(sum.fn.stk) <- recp.lat 

    } else {

        # get NEE files and fluxes based on extent of footprints
        nee.files <- list.files(nee.path, 'hrly', full.names = T)
        nee.file <- nee.files[grepl(substr(timestr, 1, 6), nee.files)]
        if (FALSE %in% file.exists(nee.file) | length(nee.file) == 0) 
            stop('NO NEE files found for this time\n')

        # read hourly NEE fluxes, given backward time
        cat('couple.xfoot.nee(): Extracting NEE fluxes based on footprints\n')
        nee.stk <- crop.smurf.hrly.flux(timestr, nee.file, site.ext = foot.ext,
                                        varname = nee.varname, nhrs = foot.nhrs, 
                                        res = foot.res)
        
        # remove the last layer (flux right at receptor time), if nhrs < 0
        nee.stk <- nee.stk[[-which(names(nee.stk) == paste0('X', timestr))]]
        #nee.stk[is.na(nee.stk)] <- 0

        # loop over each receptor 
        sum.fn.stk <- NULL
        for (f in 1 : length(foot.files)) {

            print(f)
            foot.stk <- stack(foot.files[f])
            if (nlayers(foot.stk) == 1) next 
            
            foot.names <- as.numeric(gsub('X', '', names(foot.stk)))
            foot.dates <- as.POSIXct(foot.names, origin = '1970-01-01 00:00:00',
                                     tz = 'UTC')
            foot.dates <- lubridate::floor_date(foot.dates, 'hour')
            foot.timestr <- as.numeric(paste0(format(foot.dates, '%Y%m%d%H')))
            names(foot.stk) <- foot.timestr 

            # crop footprint based on extent of nee
            # in case NEE's extent is smaller than foot's
            sel.foot <- crop(foot.stk, extent(nee.stk))

            # subset NEE fluxes based on the available hours in foots
            # sometimes # of hours of foot is smaller than that of fluxes
            sel.nee <- subset(nee.stk, which(names(nee.stk) %in% names(sel.foot)))
            foot.nee <- sel.foot * sel.nee; names(foot.nee) <- names(sel.nee)
            
            # sum the hourly bio contributions to one map 
            sum.fn <- sum(foot.nee)
            if (f == 1) sum.fn.stk <- sum.fn else sum.fn.stk <- stack(sum.fn.stk, sum.fn)
        }

        names(sum.fn.stk) <- recp.lat   # assign receptor latitude as layer name
        writeRaster(sum.fn.stk, filename = fn, overwrite = T)

    }   # end if
   
    foot.info$xco2.bio <- as.numeric(cellStats(sum.fn.stk, 'sum'))
    foot.info$timestr <- timestr

    fn <- file.path(store.path, paste0('XCO2bio_', site, '_', timestr, '.txt'))
    write.table(foot.info, file = fn, sep = ',', row.names = F, quote = F)
    cat(paste('couple.xfoot.nee(): txtfile stored as', fn, '\n\n'))

    return(list(sum.fn.stk = sum.fn.stk, foot.info = foot.info))

}   # end of subroutine


# ---------------------------------------------------------------------------- #
plot.xco2.bio <- function(site, timestr, xbio.rt, foot.info, zoom = 8, 
                          font.size = rel(0.8)) {
    
    title <- paste('Lat-int XCO2.bio [ppm] and total XCO2.bio at receptors for', 
                    site, 'on', timestr)

    loc <- suppressWarnings(get.lon.lat(site, dlon = 1, dlat = 2))
    map <- ggplot.map(map = 'ggmap', maptype = 'terrain', color = 'bw', zoom = zoom, 
                      center.lat = loc$sitelat, center.lon = loc$sitelon)
    
    b0 <- theme(legend.position = 'bottom',
                legend.text = element_text(size = font.size),
                legend.key = element_blank(), 
                legend.key.height = unit(0.5, 'cm'),
                legend.key.width = unit(1.2, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size)) 

    # convert raster to data.frame
    #xbio.rt[xbio.rt == 0] <- NA
    xbio.df  <- as.data.frame(xbio.rt, xy = T) #%>% na.omit()

    max.y1 <- max(abs(xbio.df$layer))
    b1 <- map[[1]] + coord_equal(1.2) + b0 + labs(x = 'LONGITUDE', y = 'LATITUDE') + 
          geom_raster(data = xbio.df, aes(x + map[[3]], y + map[[2]], fill = layer), alpha = 0.7) + 
          scale_fill_gradientn(name = 'Lat-integrated\ngridded XCO2.bio',
                               limits = c(-max.y1, max.y1), 
                               colors = rev(brewer.pal(11, 'BrBG'))) + 
          xlim(c(loc$minlon, loc$maxlon)) 
    
    max.y2 <- max(abs(foot.info$xco2.bio))
    b2 <- map[[1]] + b0 + labs(x = 'LONGITUDE', y = 'LATITUDE') +
          geom_point(data = foot.info, aes(lon, lat, fill = xco2.bio), 
                     shape = 21, color = 'gray70') + 
          scale_fill_gradientn(name = 'XCO2.bio\nat each X-recptor', 
                               colors = terrain.colors(20)) +
                               #rev(brewer.pal(11, 'BrBG'))) + #, 
                               #limits = c(-max.y2, max.y2)) + 
          xlim(c(loc$minlon, loc$maxlon)) 

    bb <- ggarrange(b1, b2, ncol = 2, widths = c(1.3, 1))
    bb <- annotate_figure(bb, top = title)

    return(bb)
}