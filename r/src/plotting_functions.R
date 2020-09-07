#' several plotting scripts including 
#' MODIS IGBP, GPP-SIF slopes, direct OCO-2 SIF, 
#' spatial SIF, spatial GPP, and tree density, and aboveground biomass

#' Dien Wu, 04/19/2019
#' allow for raster layer as input for plotting


def.col <- function(){
  return(c('black', '#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8',
           '#A7DA64','#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'))
}


# get coutnry and county outline for levelplot()
get.outline.levelplot <- function() {
    
    library(maptools)
    # paths for country and county outlines shapefiles
    homedir <- '/uufs/chpc.utah.edu/common/home'
    input.path <- file.path(homedir, 'lin-group7/wde/input_data')
    world.path <- file.path(input.path, 'world_shapefile/TM_WORLD_BORDERS-0.3.shp')
    us.path <- file.path(input.path, 'TIGER/tl_2018_us_county.shp')

    # read as SpatialPolygonsDataFrame
    worldmap <- readShapeSpatial(world.path)    
    uscountymap <- readShapeSpatial(us.path)

    proj4string(worldmap) <- "+proj=longlat +datum=WGS84"
    proj4string(uscountymap) <- "+proj=longlat +datum=WGS84"

    list(worldmap = worldmap, uscountymap = uscountymap)
}   # end of function


# --------------------------------------------------------------------------- #
ggmap.igbp <- function(igbp, site, yr, map, font.size = rel(0.8), res = '500m', 
                       scale.coord = 1.2, legend.ncol = 4) {

    # start plotting
    title <- paste(res, 'IGBP from MCD12Q1\nover', site, 'for year', yr)
    igbp.col.df <- get.igbp.col() # get igbp colors

    name <- igbp.col.df$name; attributes(name)$names <- igbp.col.df$val
    col <- igbp.col.df$col; attributes(col)$names <- igbp.col.df$val
    igbp.col.list <- list(name = name, val = igbp.col.df$val, col = col)

    if ('RasterLayer' %in% class(igbp)) 
        igbp.df <- raster::as.data.frame(igbp, xy = T) %>% 
                   dplyr::rename(lon = x, lat = y)  
    if ('data.frame' %in% class(igbp)) igbp.df <- igbp

    colnames(igbp.df)[3] <- 'z'
    if (is.null(map)) map <- list(ggplot(), 0, 0)

    i1 <- map[[1]] + coord_equal(scale.coord) + 
          geom_raster(data = igbp.df, aes(lon + map[[3]], lat + map[[2]], 
                                          fill = factor(z)), alpha = 0.6) + 
          scale_fill_manual(name = NULL, values = igbp.col.list$col, 
                            breaks = igbp.col.list$val, 
                            labels = igbp.col.list$name) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
          theme(legend.position = 'bottom',
                legend.text = element_text(size = font.size),
                legend.key = element_blank(), 
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(0.7, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size)) + 
          guides(fill = guide_legend(ncol = legend.ncol))

    return(i1)
}



# --------------------------------------------------------------------------- #
ggmap.lcz.type <- function(lcz.rt, site, font.size = rel(0.8), res = '120 m', 
                           scale.coord = 1.2, legend.ncol = 3) {

    # start plotting
    title <- paste(res, 'WUDAPT LCZ types over', site)
    lcz.info <- get.lcz.frac() # get igbp colors
    lcz.df <- as.data.frame(lcz.rt, xy = T)
    names(lcz.df) <- list('lon', 'lat', 'lcz.id')

    lcz.info[17, c('veg.frac.min', 'veg.frac.max')] <- NA
    lcz.df <- lcz.df %>% left_join(lcz.info, by = 'lcz.id') %>% na.omit()

    # create color list for l1
    name <- lcz.info$lcz.name; attributes(name)$names <- lcz.info$lcz.id
    col <- lcz.info$lcz.col; attributes(col)$names <- lcz.info$lcz.id
    lcz.col.list <- list(name = name, val = lcz.info$lcz.id, col = col)
    
    l1 <- ggplot() + coord_equal(scale.coord) + theme_bw() +
          geom_raster(data = lcz.df, aes(lon, lat, fill = factor(lcz.id)), alpha = 0.8) + 
          scale_fill_manual(name = NULL, values = lcz.col.list$col, 
                            breaks = lcz.col.list$val, labels = lcz.col.list$name) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
          theme(legend.position = 'bottom',
                legend.text = element_text(size = font.size),
                legend.key = element_blank(), 
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(0.7, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size)) + 
          guides(fill = guide_legend(ncol = legend.ncol))


    # start plotting
    title <- paste('Minimum impervious and vegetated fractions [%] for', 
                   site, 'according to WUDAPT LCZ')

    f1 <- ggplot() + coord_equal(scale.coord) + theme_bw() + 
          geom_raster(data = lcz.df, aes(lon, lat, fill = imp.frac.min), alpha = 0.8) + 
          scale_fill_gradientn(name = 'Impervious\nfraction [%]', colors = brewer.pal(9, 'YlOrRd')) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE') + 
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(0.7, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))
    
    f2 <- ggplot() + coord_equal(scale.coord) + theme_bw() + 
          geom_raster(data = lcz.df, aes(lon, lat, fill = veg.frac.min), alpha = 0.8) + 
          scale_fill_gradientn(name = 'Vegetated\nfraction [%]', colors = brewer.pal(9, 'BuGn')) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE') + 
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(0.7, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))

    ff <- ggarrange(f1, f2, ncol = 2)
    ff <- annotate_figure(ff, top = text_grob(title))
    lff <- ggarrange(l1, ff, ncol = 2, widths = c(1, 1.7))
    return(lff)
}



# --------------------------------------------------------------------------- #
ggmap.lcz.frac <- function(frac, site, map, 
                           frac.name = c('Vegetated', 'Impervious', 'Difference in')[1], 
                           frac.res = '120 m', frac.max = 100, frac.source = 'WUDAPT',
                           font.size = rel(0.8), scale.coord = 1.2) {

    if ('RasterLayer' %in% class(frac)) 
        frac.df <- as.data.frame(frac, xy = T) %>% dplyr::rename(lon = x, lat = y)  
    if ('data.frame' %in% class(frac)) frac.df <- frac
    colnames(frac.df)[3] <- 'z'

    frac.min <- 0
    if (frac.name == 'Vegetated') {
        col <- 'BuGn'
    } else if (frac.name == 'Impervious') {
        col <- 'YlOrRd'
    } else if (grepl('Difference', frac.name)) {
        col <- 'BrBG'
        frac.min <- -frac.max
    }

    f1 <- map[[1]] + coord_equal(scale.coord) + 
          geom_raster(data = frac.df %>% na.omit(), 
                      aes(lon + map[[3]], lat + map[[2]], fill = z), alpha = 0.7) + 
          scale_fill_gradientn(name = paste(frac.name, '\nfraction [%]'), 
                               colors = brewer.pal(9, col), 
                               limit = c(frac.min, frac.max)) + 
          #xlim(c(min(frac.df$lon), max(frac.df$lon))) + 
          #ylim(c(min(frac.df$lat), max(frac.df$lat))) +
          
          labs(x = 'LONGITUDE', y = 'LATITUDE', 
               title = paste(frac.res, frac.name, 'fraction [%] from', 
                             frac.source)) + #, 'for', site)) + 
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(0.7, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))

    return(f1)
}



# --------------------------------------------------------------------------- #
# slope can only be in form of data frame!!!
ggmap.gpp.sif.slope <- function(slope, site, yr, map, font.size = rel(0.8), 
                                ref = c('Sun et al., 2018', 'Zhang et al., 2018')[2]) {

    # plot GPP-SIF slopes
    title <- paste('500m empirical GPP-SIF slopes [gC/m2/d GPP : W/m2/um/sr SIF]',
                   '\nfrom', ref, 'over', site, 'for year', yr)

    if ('RasterLayer' %in% class(slope)) 
        stop('Incorrect form of initial slopes, cannot be raster; need IGBP')
    if ('data.frame' %in% class(slope)) slope.df <- slope

    s1 <- map[[1]] + coord_equal(1.2) + 
          geom_raster(data = slope.df %>% na.omit(), 
                      aes(lon + map[[3]], lat + map[[2]], 
                          fill = factor(paste(SLP, abbr))), alpha = 0.4) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title, 
               fill = 'GPP-SIF slopes\nfor IGBP') + 
          theme(legend.position = 'bottom',
                legend.text = element_text(size = font.size),
                legend.key = element_blank(), 
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(1, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size)) + 
        guides(fill = guide_legend(ncol = 3))
    return(s1)
}


# --------------------------------------------------------------------------- #
ggmap.gpp.sif.slope.coarse <- function(slope, site, timestr, map,
                                       font.size = rel(1.0), res = '0.05◦', 
                                       zlim = NULL){

    # plot 0.05 deg slopes on map
    if ('RasterLayer' %in% class(slope)) 
        slope.df <- raster::as.data.frame(slope, xy = T) %>% 
                    rename(lon = x, lat = y)
    if ('data.frame' %in% class(slope)) slope.df <- slope
    colnames(slope.df)[3] <- 'layer'

    if (is.null(zlim)) zlim <- c(min(slope.df$layer), max(slope.df$layer))
    title <- paste(res, 'GPP-SIF slopes over', site, 'for', timestr)
    s2 <- map[[1]] + coord_equal(1.2) + 
          geom_raster(data = slope.df %>% na.omit(), 
                      aes(lon + map[[3]], lat + map[[2]], fill = layer), alpha = 0.4) + 
          scale_fill_gradientn(name = 'GPP:SIF slopes', colours = ggdef.col(n = 10), 
                               limits = zlim) +
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(1, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))
    return(s2)
}


# --------------------------------------------------------------------------- #
# map can be generated from calling ggplot.map()
ggmap.sif <- function(oco2.sif, site, timestr, map, font.size = rel(0.8)){

    title <- paste('Direct OCO-2 SIF > -0.2 [W/m2/sr/µm]\naveraged over 757 and 771nm for', 
                    site, 'on', timestr)
    o1 <- map[[1]] + coord_equal(1.2) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) +
          geom_point(data = oco2.sif %>% filter(avg.sif > -0.2) %>% na.omit(), 
                     aes(lon + map[[3]], lat + map[[2]], colour = avg.sif), size = 0.2) +
          scale_colour_gradientn(name = 'OCO-2 SIF', colours = rev(terrain.colors(10))) +
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(1, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))
    return(o1)
}


# --------------------------------------------------------------------------- #
# map can be generated from calling ggplot.map()
ggmap.csif <- function(csif, site, timestr, map = NULL, res = '0.05◦', 
                       font.size = rel(0.8), zlim = NULL, var.title = 'CSIF', 
                       scale.coord = 1.2){

    # convert rasterlayer to data.frame if needed
    if ('RasterLayer' %in% class(csif)) 
        csif.df <- raster::as.data.frame(csif, xy = T) %>% 
                   dplyr::rename(lon = x, lat = y)

    if ('data.frame' %in% class(csif)) csif.df <- csif

    # plot CSIF on map and overlay with direct SIF measurements from OCO2
    ini.z <- colnames(csif.df)[3]; colnames(csif.df)[3] <- 'z'

    if (is.null(zlim)) zlim <- c(min(csif.df$z), max(csif.df$z))
    title <- paste(res, var.title, 'over', site, 'on', timestr)

    if (is.null(map)) {
        map <- list(ggplot() + theme_bw(), 0, 0)
    } else if ('list' %in% class(map)) {
        map[[1]] <- map[[1]] + coord_equal(scale.coord)
    } else if (!'list' %in% class(map)) map <- list(map, 0, 0)

    c1 <- map[[1]] + geom_raster(data = csif.df %>% na.omit(), 
                     aes(lon + map[[3]], lat + map[[2]], fill = z), alpha = 0.9) + 
          scale_fill_gradientn(name = ini.z, limits = zlim, 
                               colours = rev(terrain.colors(10))) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(1, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))
    return(c1)
}


# --------------------------------------------------------------------------- #
# map can be generated from calling ggplot.map()
ggmap.gosif <- function(gosif, site, timestr, map, font.size = rel(0.8), zlim = NULL){

    # convert rasterlayer to data.frame if needed
    if ('RasterLayer' %in% class(gosif))
        gosif.df <- raster::as.data.frame(gosif, xy = T) %>% 
                    rename(lon = x, lat = y, GOSIF = layer)
    if ('data.frame' %in% class(gosif)) gosif.df <- gosif

    # plot GOSIF on map and overlay with direct SIF measurements from OCO2
    if (is.null(zlim)) zlim <- c(min(gosif.df$GOSIF), max(gosif.df$GOSIF))

    title <- paste('0.05◦ 8-day GOSIF [Li and Xiao., 2019]\nover', site, 'for', timestr)
    g1 <- map[[1]] + coord_equal(1.2) + 
            geom_raster(data = gosif.df %>% na.omit(), 
                        aes(lon + map[[3]], lat + map[[2]], fill = GOSIF), alpha = 0.4) + 
            scale_fill_gradientn(name = 'GOSIF', limits = zlim, 
                                 colours = rev(terrain.colors(10))) + 
            labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
            theme(legend.position = 'bottom', legend.key = element_blank(), 
                  legend.text = element_text(size = font.size),
                  legend.key.height = unit(0.3, 'cm'),
                  legend.key.width = unit(1, 'cm'),
                  axis.title.y = element_text(size = font.size, angle = 90),
                  axis.title.x = element_text(size = font.size, angle = 0),
                  axis.text = element_text(size = font.size),
                  axis.ticks = element_line(size = font.size),
                  title = element_text(size = font.size),
                  strip.text = element_text(size = font.size))
    return(g1)
}


# --------------------------------------------------------------------------- #
# if res is coarse (only plotting few grids), please set `map` to NULL
# otherwise plotting coarse grid on map will make the map distorted
ggmap.gpp <- function(gpp, site, timestr, map = NULL, res = '0.05◦', 
                      zlim = NULL, font.size = rel(0.8), unit = '[umol/m2/s]', 
                      type = c('CSIF_clear', 'CSIF_all', 'GOSIF', 'MsTMIP')[1], 
                      zname = 'GPP', scale.coord = 1.2) {

    # convert rasterlayer to data.frame if needed
    if ('RasterLayer' %in% class(gpp)) 
        gpp.df <- raster::as.data.frame(gpp, xy = T) %>% 
                  dplyr::rename(lon = x, lat = y) %>% na.omit()
    if ('data.frame' %in% class(gpp)) gpp.df <- gpp %>% na.omit()

    colnames(gpp.df)[3] <- 'layer'
    if (is.null(zlim)) zlim <- c(min(gpp.df$layer), max(gpp.df$layer))
    title <- paste(res, 'GPP', unit, 'based on', type, '\nover', site, 'on', timestr)

    if (is.null(map)) {
        map <- list(ggplot() + theme_bw(), 0, 0)
    } else if ('list' %in% class(map)) {
        map[[1]] <- map[[1]] + coord_equal(scale.coord)
    } else if (!'list' %in% class(map)) map <- list(map, 0, 0)

    g1 <- map[[1]] + geom_raster(data = gpp.df, aes(lon + map[[3]], lat + map[[2]], 
                                                    fill = layer), alpha = 0.6) +
          scale_fill_gradientn(name = zname, limits = zlim, 
                               colours = rev(terrain.colors(10))) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.6, 'cm'),
                legend.key.width = unit(1.5, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))

    return(g1)
}


# --------------------------------------------------------------------------- #
# plot hourly GPP, must be RasterLayer
ggmap.hrly.gpp <- function(gpp.rt, site, timestr, map, 
                           type = c('CSIF_clear', 'CSIF_all', 'GOSIF', 'MsTMIP')[1],
                           res = 0.05, unit = '[umol/m2/s]') {

    gppTheme <- rasterTheme(region = rev(terrain.colors(12)))
    min.gpp <- min(getValues(gpp.rt), na.rm = T)
    max.gpp <- max(getValues(gpp.rt), na.rm = T)
    gpp.break <- max.gpp / 15 # 15 breaks
    gppLab <- seq(0, max.gpp, gpp.break)

    title <- paste(res, '◦ hourly GPP', unit, 'based on', type, 
                        'over', site, 'for', timestr)
    g1 <- levelplot(gpp.rt, par.settings = gppTheme, margin = F, at = gppLab,
                    alpha.regions = 0.7, colorkey = list(space = 'right'), 
                    main = title, maxpixels = 1E6, 
                    xlab = 'LONGITUDE', ylab = 'LATITUDE') 

    return(g1)
}


# --------------------------- tree density ---------------------------------- #
# map can be generated from calling ggplot.map()
ggmap.tree.density <- function(td.df, site, timestr, map, font.size = rel(1.0), 
                               res = '1/120') {

    # plot 0.05 deg GPP on map
    title <- paste(res, '◦ Tree Density over', site, 'for', timestr)
    t1 <- map[[1]] + coord_equal(1.2) + 
          geom_raster(data = td.df %>% na.omit(), 
                      aes(lon + map[[3]], lat + map[[2]], fill = layer), alpha = 0.4) + 
          scale_fill_gradientn(name = 'Tree Density [#/grid]', trans = 'log10',
                               colours = rev(terrain.colors(10))) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(1, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))
    return(t1)
}


# ------------------------- AG biomass data --------------------------------- #
# map can be generated from calling ggplot.map()
ggmap.agb <- function(agb, site, map, font.size = rel(0.8), 
                      agb.type = c('DAAC_ORNL', 'GEOCARBON', 'GlobBiomass')[3], 
                      agb.res = c('1◦', '0.01◦', '~100 m')[3], 
                      agb.unit = 'tonnes/ha', scale.coord = 1.2) {
    
    # convert rasterlayer to data.frame if needed
    if ('RasterLayer' %in% class(agb)) 
        agb.df <- raster::as.data.frame(agb, xy = T) %>% 
                  dplyr::rename(lon = x, lat = y) %>% na.omit()
    if ('data.frame' %in% class(agb)) agb.df <- agb %>% na.omit()
    colnames(agb.df)[3] <- 'AGB'

    a1 <- map[[1]] + coord_fixed(scale.coord) +
          geom_raster(data = agb.df, aes(lon + map[[3]], lat + map[[2]], 
                      fill = AGB), alpha = 0.5) + 
          scale_fill_gradientn(colours = rev(terrain.colors(10)), 
                               name = paste0('AGB\n[', agb.unit, ']')) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', 
               title = paste(agb.res, agb.type, 'AGB for', site)) + 
          theme(legend.position = 'bottom', legend.key = element_blank(),
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(1, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                    axis.title.x = element_text(size = font.size, angle = 0),
                    axis.text = element_text(size = font.size),
                    axis.ticks = element_line(size = font.size),
                    title = element_text(size = font.size),
                    strip.text = element_text(size = font.size))
    return(a1)
}


# -------------------------------------------------------------------------
# plot hours as centered hrs, 0.5 for 00-01UTC
ggplot.sf.rad <- function(sf.cru.df, sf.clm.df, timestr, site, font.size = rel(1.0)) {

    t1 <- ggplot() + theme_bw() + 
          labs(x = 'Centered UTC [hr]', y = 'SW radiation scaling factors', 
               title = paste('Time series of scaling factors from CRU (6hrly, 0.5x0.5deg)',
                             '\nwith CLM output (hrly, 0.625x0.47deg)\nfor', site, 
                             'on', timestr)) +
          geom_point(data = sf.clm.df, aes(hr + 0.5, sf, colour = 'CLM'), size = 3) + 
          geom_line(data = sf.clm.df, aes(hr + 0.5, sf, colour = 'CLM'), linetype = 2) + 
          geom_point(data = sf.cru.df, aes(hr + 3, sf, colour = 'CRU'), 
                     shape = 17, size = 3) + 
          scale_x_continuous(breaks = seq(0, 24, 2), labels = seq(0, 24, 2)) + 

          theme(legend.position = 'bottom', legend.key = element_blank(),
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(1, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                    axis.title.x = element_text(size = font.size, angle = 0),
                    axis.text = element_text(size = font.size),
                    axis.ticks = element_line(size = font.size),
                    title = element_text(size = font.size),
                    strip.text = element_text(size = font.size))

    return(t1)
}

# -------------------------------------------------------------------------
# plot difference in gridded data 
ggmap.diff <- function(diff, site, timestr, map = NULL, text, 
                       font.size = rel(0.8), scale.coord = 1.2) {

    # convert rasterlayer to data.frame if needed
    if ('RasterLayer' %in% class(diff)) 
        diff.df <- raster::as.data.frame(diff, xy = T) %>% 
                   dplyr::rename(lon = x, lat = y) %>% na.omit()
    if ('data.frame' %in% class(diff)) diff.df <- diff %>% na.omit()
    colnames(diff.df)[3] <- 'layer'

    zlim <- c(-max(abs(diff.df$layer), na.rm = T), 
               max(abs(diff.df$layer), na.rm = T))  # zero -- white color

    title <- paste('Difference in', text, '\nover', site, 'for', timestr)

    if (is.null(map)) {
        map <- list(ggplot() + theme_bw(), 0, 0)
    } else if ('list' %in% class(map)) {
        map[[1]] <- map[[1]] + coord_equal(scale.coord)
    } else if (!'list' %in% class(map)) map <- list(map, 0, 0)

    d1 <- map[[1]] + geom_raster(data = diff.df, aes(lon, lat, fill = layer), 
                                 alpha = 0.9) + 
          scale_fill_gradientn(name = 'DIFF', colours = rev(brewer.pal(11, 'RdBu')), 
                               limits = zlim) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
             
          theme(legend.position = 'bottom', legend.key.width = unit(1, 'cm'),
                legend.text = element_text(size = font.size),
                legend.key = element_blank(), legend.key.height = unit(0.3, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))
    return(d1)
}


# -------------------------------------------------------------------------
# plot tower Reco with potential features 
#all.f <- c('GPP', 'Tair', 'Tsoil', 'rH', 'VPD', 'DoY')
plot.reco.feature.tower <- function(obs.df, 
                                    features = c('GPP', 'Tair', 'Tsoil', 'rH', 
                                                 'VPD', 'DoY'), 
                                    y = 'Reco', colour = 'DoY', title = NULL, 
                                    nrow = 2, ncol = 3, 
                                    font.size = rel(0.7)) {

    # plot daily mean Reco with all features before training
    sz <- 0.05; col <- c(brewer.pal(9, 'YlOrBr'), rev(brewer.pal(9, 'YlOrBr'))[-1])

    p0 <- ggplot(obs.df, aes_string(y = y, colour = colour)) + 
          theme(legend.position = 'bottom', legend.key = element_blank(),
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.2, 'cm'), legend.key.width = unit(3, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size)) + 
            scale_colour_gradientn(colours = col, breaks = seq(30, 360, 90), 
                                   labels = seq(30, 360, 90)) + theme_classic() 
    if (!is.null(title)) p0 <- p0 + labs(title = title)
    
    plot.list <- list()
    for (f in 1 : length(features)) {
        xseq <- 5; yseq <- 2; lim1 <- -40; lim2 <- 100
        if (grepl('norm', features[f])) { xseq <- 0.2; yseq <- 2; lim1 <- 0; lim2 <- 1 }
        if (grepl('cumsum', features[f])) { xseq <- 20; lim2 <- 500 }
        
        p1 <- p0 + geom_point(aes_string(x = features[f]), size = sz) + 
                   scale_x_continuous(breaks = seq(lim1, lim2, xseq), 
                                      labels = seq(lim1, lim2, xseq)) + 
                   scale_y_continuous(breaks = seq(lim1, lim2, yseq), 
                                      labels = seq(lim1, lim2, yseq)) 

        if (grepl('GPP', features[f]) & !grepl('cumsum', features[f]))
            p1 <- p1 + geom_line(aes_string(x = y), linetype = 2, colour = 'gray50') 
                                
        plot.list[[f]] <- p1
    }
        
    panel.lab <- c('a)', 'b)', 'c)', 'e)', 'f)', 'g)', 'i)', 'j)', 'k)')[1: length(features)]
    pp <- ggarrange(plotlist = plot.list, nrow = nrow, ncol = ncol, common.legend = T,#) 
                    labels = panel.lab)
    return(pp)
}


# -------------------------------------------------------------------------
# plot tower GPP with potential features 
#all.f <- c('GPP', 'Tair', 'Tsoil', 'rH', 'VPD', 'DoY')
plot.gpp.feature.tower <- function(obs.df, font.size = rel(0.7), 
                                   features = c('SIF', 'Tair', 'DoY'), 
                                   y = 'GPP', colour = 'DoY', nrow = 1, ncol = 3) {

    # plot daily mean Reco with all features before training
    sz <- 0.05; col <- c(brewer.pal(9, 'Oranges'), rev(brewer.pal(9, 'Oranges'))[-1])

    p0  <-  ggplot(obs.df, aes_string(y = y, colour = colour)) + 
            theme(legend.position = 'bottom', legend.key = element_blank(),
                  legend.text = element_text(size = font.size),
                  legend.key.height = unit(0.3, 'cm'), 
                  legend.key.width = unit(1.5, 'cm'),
                  axis.title.y = element_text(size = font.size, angle = 90),
                  axis.title.x = element_text(size = font.size, angle = 0),
                  axis.text = element_text(size = font.size),
                  axis.ticks = element_line(size = font.size),
                  title = element_text(size = font.size),
                  strip.text = element_text(size = font.size)) + 
            theme_classic() + scale_colour_gradientn(colours = col) 

    plot.list <- list()
    for (f in 1 : length(features)) 
        plot.list[[f]] <- p0 + geom_point(aes_string(x = features[f]), size = sz) + 
                               geom_smooth(aes_string(x = features[f]), 
                                           method = 'lm', se = FALSE)

    pp <- ggarrange(plotlist = plot.list, nrow = nrow, ncol = ncol, 
                    common.legend = T)
    return(pp)
}



# --------------------------------------------------------------------------- #
# map can be generated from calling ggplot.map()
# plotting Tair or Tsoil
ggmap.temp <- function(temp, site, timestr, map, font.size = rel(1.0), 
                       zlim = NULL, param = c('Tair', 'Tsoil'), by = 5 ){

    # convert rasterlayer to data.frame if needed
    if ('RasterLayer' %in% class(temp))
        temp.df <- raster::as.data.frame(temp, xy = T) %>% rename(lon = x, lat = y)
    if ('data.frame' %in% class(temp)) temp.df <- temp

    # plot gridded temp on map and overlay with direct SIF measurements from OCO2
    colnames(temp.df)[3] <- 'z'

    if (is.null(zlim)) zlim <- c(min(temp.df$z), max(temp.df$z))
    title <- paste(param, '◦C over', site, 'for', timestr)

    t1 <- map[[1]] + coord_equal(1.2) + 
          geom_raster(data = temp.df %>% na.omit(), 
                    aes(lon + map[[3]], lat + map[[2]], fill = z), alpha = 0.4) + 
          scale_fill_gradientn(name = 'Temperature [◦C]', limits = zlim, 
                               colours = rev(heat.colors(10)), 
                               breaks = seq(-40, 50, by), labels = seq(-40, 50, by)) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.3, 'cm'),
                legend.key.width = unit(1, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))
    return(t1)
}



# --------------------------------------------------------------------------- #
ggmap.reco <- function(reco, site, timestr, map, res = '0.05◦', zlim = NULL, 
                       font.size = rel(0.8), reco.unit = '[umol/m2/s]', 
                       scale.coord = 1.2) {

    # convert rasterlayer to data.frame if needed
    if ('RasterLayer' %in% class(reco) | 'RasterBrick' %in% class(reco)) 
        reco.df <- raster::as.data.frame(reco, xy = T) %>% 
                   dplyr::rename(lon = x, lat = y) %>% na.omit()
    if ('data.frame' %in% class(reco)) reco.df <- reco %>% na.omit()

    colnames(reco.df)[3] <- 'layer'
    if (is.null(zlim)) zlim <- c(min(reco.df$layer), max(reco.df$layer))
    title <- paste(res, 'Reco', reco.unit, 'over', site, '\non', timestr)

    r1 <- map[[1]] + coord_equal(scale.coord) + 
          geom_raster(data = reco.df, aes(lon + map[[3]], lat + map[[2]], 
                                          fill = layer), alpha = 0.6) 

    r1 <- r1 + scale_fill_gradientn(name = 'Reco', limits = zlim,
                                    colours = rev(terrain.colors(10))) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.6, 'cm'),
                legend.key.width = unit(1.5, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))

    return(r1)
}



# --------------------------------------------------------------------------- #
ggmap.nee <- function(nee, site, timestr, map, res = '0.05◦', zlim = NULL, 
                       font.size = rel(0.8), nee.unit = '[umol/m2/s]') {

    # convert rasterlayer to data.frame if needed
    if ('RasterLayer' %in% class(nee)) nee.df <- raster::as.data.frame(nee, xy = T) %>% 
                                                 rename(lon = x, lat = y) %>% na.omit()
    if ('data.frame' %in% class(nee)) nee.df <- nee %>% na.omit()

    colnames(nee.df)[3] <- 'layer'
    if (is.null(zlim)) zlim <- c(-max(abs(nee.df$layer)), max(abs(nee.df$layer)))
    title <- paste(res, 'NEE', nee.unit, 'over', site, '\non', timestr)

    n1 <- map[[1]] + coord_equal(1.2) + 
          geom_raster(data = nee.df, aes(lon + map[[3]], lat + map[[2]], 
                                         fill = layer), alpha = 0.6) 

    n1 <- n1 + scale_fill_gradientn(name = 'NEE', limits = zlim,
                                    colours = rev(brewer.pal(11, 'RdBu'))) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
          theme(legend.position = 'bottom', legend.key = element_blank(), 
                legend.text = element_text(size = font.size),
                legend.key.height = unit(0.6, 'cm'),
                legend.key.width = unit(1.5, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))

    return(n1)
}

# end of scripts



convert.loc2spdf <- function(sites, dlat = 2, dlon = 2){

    loc <- suppressWarnings(get.lon.lat(sites, dlat = dlat, dlon = dlon))
    ll.ref <- loc %>% dplyr::select(site, minlon, maxlon, minlat, maxlat)
    ll.rev <- NULL
    for (s in 1:length(sites)) {
        tmp <- ll.ref %>% filter(site == sites[s])
        tmp.ll  <- data.frame(lon = c(tmp$minlon, tmp$minlon, tmp$maxlon, tmp$maxlon),
                              lat = c(tmp$minlat, tmp$maxlat, tmp$maxlat, tmp$minlat),
                              site = rep(sites[s], 4))
        ll.rev <- rbind(ll.rev, tmp.ll)
    }

    # make a list, and only select lon, lat,
    # ref from https://gis.stackexchange.com/questions/171124/data-frame-to-spatialpolygonsdataf
    # rame-with-multiple-polygons
    ll.list <- split(ll.rev, ll.rev$site)
    ll.list <- lapply(ll.list, function(x) { x['site'] <- NULL; x })

    # convert to polygon and add id variable
    ps <- lapply(ll.list, Polygon)
    ps <- lapply(seq_along(ps), function(i) Polygons(list(ps[[i]]),
                                            ID = names(ll.list)[i]))

    # create SpatialPolygons object
    ll.sp <- SpatialPolygons(ps, proj4string = CRS('+proj=longlat +datum=WGS84') )

    # convert to SP dataframe
    centroids <- coordinates(ll.sp)
    x <- centroids[,1]
    y <- centroids[,2]
    SpatialPolygonsDataFrame(ll.sp, data = data.frame(x, y, row.names = row.names(ll.sp)))
}
