#' all get.*.r functions, DW, 08/05/2019
#' @author Dien Wu
#' get_* functions including get.lon.lat(), get.igbp.abbr(), get.igbp.col(), 
#'                           get.raw.slp(), get.fluxnet.info(), get.reco.cv(), 


# ---------------------------------------------------------------------------- #
# functions to create the latitude dependence of tree type fractions
# DW, 03/28/2020 
# ---------------------------------------------------------------------------- #
# return fractions of DBF, EBF, and ENF at a given latitude [-60, 70] in a data frame
get.tree.type.frac <- function(urb.lat, smurf_wd = smurf_wd) {

    # get zonal mean relative tree fractions from txt file 
    zonal.file <- file.path(smurf_wd, 'data/zonal_tree_frac_globe.txt')
    zonal.df <- read.table(zonal.file, sep = ',', header = T)

    uni.lat <- unique(zonal.df$lat) # unique latitude in zonal mean extraction
    match.lat <- uni.lat[findInterval(urb.lat, uni.lat)]    # match lats

    frac.df <- data.frame(urb.lat = urb.lat, match.lat = match.lat) %>% 
              left_join(zonal.df, by = c('match.lat' = 'lat'))

    # sanity check 
    if (F) {
        frac.long <- melt(frac.df %>% dplyr::select(-c(match.lat)), id.vars = 'urb.lat')
        f1 <- ggplot(data = frac.long) + theme_bw() + 
              geom_point(aes(urb.lat, value, color = variable)) 
            
    }   # end if

    return(frac.df) 
}   # end of get.tree.type.frac


# ---------------------------------------------------------------------------- #
#' script to define spatial domain
#' @author Dien Wu, 07/05/2018
# ---------------------------------------------------------------------------- #
#' update:
#' use geocode and SpatialPoints to find lat/lon coordinates,
#' country and reg name, DW, DR, 08/15/2018
#' site can be a vector
#' add coordinates for flux tower, rename output variables, DW, 05/21/2018
#' flux.file downloaded from FLUXNET 2015 
#' add 'trunc' to lat/lon domain, DW, 06/10/2019 

get.lon.lat <- function(site, dlon, dlat, loc = NULL, 
    flux.file = file.path('/uufs/chpc.utah.edu/common/home/lin-group7/wde/input_data',
                          'FLUXNET/data/fluxnet_site_info_all.csv')) {

    library(ggmap); library(rworldmap); library(sp); library(lutz)

    # if `- `is contained in the site name, site is a flux tower
    fluxTF <- grepl('-', site)   

    if (TRUE %in% fluxTF) { # if site is a FLUXNET site name
        # get site info, e.g., lat, lon 
        flux.info <- read.csv(flux.file, header = T, sep = ',', stringsAsFactors = F) %>% 
                    filter(fluxnetid %in% site)
        loc <- data.frame(lon = flux.info$longitude, lat = flux.info$latitude)
    } # end if fluxTF

    # location name to lon, lat coordinates
    if (is.null(loc)) loc <- geocode(location = site, output = 'latlon', 
                                    source = 'google', override_limit = T)

    # from https://stackoverflow.com/questions/21708488/
    # get-country-and-continent-from-longitude-and-latitude-point-in-r
    # use high res map from rworldxtra if you were concerned about detail
    countriesSP <- getMap(resolution = 'low')

    # converting points to a SpatialPoints object
    # setting CRS directly to that from rworldmap
    pointsSP <- SpatialPoints(loc, proj4string = CRS(proj4string(countriesSP)))

    # use 'over' to get indices of the Polygons object containing each point
    indices <- over(pointsSP, countriesSP)

    # get time.zone 
    tz <- tz_lookup_coords(loc$lat, loc$lon)

    # convert indices from factors to characters
    lon.lat <- data.frame(site = site, sitelon = loc$lon, sitelat = loc$lat, tz, 
                            countryid = as.character(indices$ADMIN), 
                            regid = as.character(indices$continent), 
                            iso3 = as.character(indices$ISO3),
                            minlon = loc$lon - dlon, 
                            maxlon = loc$lon + dlon,
                            minlat = loc$lat - dlat, 
                            maxlat = loc$lat + dlat, 
                            stringsAsFactors = F)

    return(lon.lat)
}   # end of get.lon.lat()


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
get.igbp.abbr <- function(input.full.name = NULL) {

    full.names <- c('Cropland', 'Croplands', 'Cropland/Natural Vegetation Mosaic', 
                    'Closed Shrublands', 
                    'Deciduous Broadleaf Forest', 'Deciduous Broadleaf Forests', 
                    'Deciduous Needleleaf Forest', 'Deciduous Needleleaf Forests',
                    'Evergreen Broadleaf Forest', 'Evergreen Broadleaf Forests',
                    'Evergreen Needleleaf Forest', 'Evergreen Needleleaf Forests', 
                    'Grassland', 'Grasslands', 
                    'Mixed Forest', 'Mixed Forests', 
                    'Open Shrublands', 
                    'Savanna', 'Savannas', 
                    'Urban and Built-Up', 
                    'Wetlands', 'Permanent Wetlands',
                    'Woody Savannas', 'Snow and Ice', 'Water')

    abbrs <- c(rep('CRO', 3), 'CSHR', rep('DBF', 2), rep('DNF', 2), 
               rep('EBF', 2), rep('ENF', 2), rep('GRA', 2), rep('MF', 2), 'OSHR', 
               rep('SAV', 2), 'URB', rep('WET', 2), 'WSAV', 'SI', 'WAT')

    df <- data.frame(name = full.names, abbr = abbrs, stringsAsFactors = F)

    if ( is.null(input.full.name)) return(df)
    if (!is.null(input.full.name)) {
        find.abbr <- df$abbr[match(input.full.name, df$name)]
        return(find.abbr)
    }   # end if 

}   # end of function


# ---------------------------------------------------------------------------- #
# get IGBP colors first ,according to MCD12 product
# add C3 /C4 CRO as additional categories
# ---------------------------------------------------------------------------- #
get.igbp.col <- function() {

    val  <- c(seq(1, 17), 254, 20, 21, 22)
    name <- c('ENF', 'EBF',   # Evergreen Needleleaf/Broadleaf Forest
              'DNF', 'DBF',   # Deciduous Needleleaf/Broadleaf Forest
              'MF',           # Mixed forest
              'CSHR', 'OSHR', # Open/Closed shrubland
              'WSAV', 'SAV', 'GRA',  # (woody) Savannas, Grasslands
              'WET',  'CRO', 'URB',  # wetland, Croplands
              'CRO/NVM', 'SI',       # Natural Vegetation Mosaic; snow and ice
              'BAR',  'WAT', 'UNC',  # barren; water; unclassified
              'WGT', 'CRO_C3', 'CRO_C4')  
              # weighted -- additional types for SMUrF ("urban trees")
              # add C3 and C4 crops

    col <- c('#008000', '#00FF00', '#99CC00', '#99FF99', '#339966', 
             '#993366', '#FFCC99', '#CCFFCC', '#FFCC00', '#FF9900', 
             '#006699', '#FFFF00', '#FF0000', '#999966', '#FFFFFF', 
             '#808080', '#000080', '#000000', '#4DC04D', '#FFFF00', '#FFFF00')
    
    out <- data.frame(val = val, name = name, col = col, stringsAsFactors = F)
    #out <- list(val = val, name = name, col = col)
    return(out)
}



# ---------------------------------------------------------------------------- #
# subroutines to assign slopes to each IGBP type 
# by Dien Wu
# ---------------------------------------------------------------------------- #
# DW, 03/29/2020: partition between C3 and C4 crops

get.raw.slp <- function(slp.file = file.path(smurf_wd, 'data/GPP-CSIF_stat_Wu.txt')) {

    # mean slopes and uncertainties for clear-sky CSIF in Zhang et al., 2018
    # for cropland, we get the slopes by training FLUXNET GPP vs. CSIF
    # treat CRO/NVM as CRO
    slope.df <- read.table(slp.file, sep = ',', header = T, stringsAsFactors = F) 
    
    # additional land types
    # WGT stands for weighted urban tree type
    # pred.slpv2() will calculate the weighted slopes, int, cv for those urban trees
    add.igbp <- c('URB', 'SI', 'BAR', 'WAT', 'UNC', 'WGT')
    add.slp <- add.int <- add.cv <- rep(0, length(add.igbp))
    add.slope.df <- data.frame(igbp = add.igbp, mean.slp = add.slp, 
                               mean.int = add.int, cv = add.cv, 
                               gpp.unit = unique(slope.df$gpp.unit), 
                               csif.unit = unique(slope.df$csif.unit), 
                               stringsAsFactors = F)
    slope.df <- rbind(slope.df, add.slope.df)

    # change CRO/NVM to CRO, add CRO_C3 and CRO_C4
    igbp.df <- get.igbp.col() %>% 
               mutate(abbr = ifelse(name == 'CRO/NVM', 'CRO', name)) %>% 
               left_join(slope.df, by = c('abbr' = 'igbp')) %>% 
               dplyr::rename(SLP = mean.slp, INT = mean.int)
    
    return(igbp.df)
}   # end of get.raw.slp()



# ---------------------------------------------------------------------------- #
# get Local Climate Zone info from WUDAPT
# ---------------------------------------------------------------------------- #
get.lcz.frac <- function() {
    lcz.id <- c(seq(1, 10, 1), seq(101, 107)) 
    lcz.name <- c('compact_high-rise', 'compact_midrise', 'compact_low-rise', 
                  'open_high-rise', 'open_midrise', 'open_low-rise', 
                  'lightweight_low-rise', 'large_low-rise', 'sparsely_built', 
                  'heavy_industry', 'dense_trees', 'scattered_trees', 'bush_scrub',
                  'low_plants', 'bare_rock_paved', 'bare_soil_sand', 'water')
    imp.frac.min <- c(40, 40, 40, rep(20, 3), 60, 30, 10, 20, rep(0, 7))
    imp.frac.max <- c(60, 70, 70, rep(40, 3), 90, 50, 20, 30, rep(10, 7))
    veg.frac.max <- c(10, 20, 30, 40, 40, 60, 30, 20, 80, 50, rep(100, 4), 10, 100, 100)
    veg.frac.min <- c(0, 0, 0, 30, 20, 30, 0, 0, 60, 40, 90, 90, 90, 90, 0, 90, 90)
    tree.frac.min <- c(rep(0, 10), 90, 90, rep(0, 5))
    tree.frac.max <- c(rep(0, 10), 100, 100, rep(0, 5))

    lcz.col <- c('#8c0000', '#d10000', '#ff0100', '#be4d01', '#ff6602', '#ff9955', 
                 '#faee05', '#bcbcbc', '#ffccaa', '#555555', '#006a01', '#01aa00', 
                 '#648526', '#b9db79', '#000000', '#fbf7ae', '#6a6aff')
    lcz.df <- data.frame(lcz.id, lcz.name, lcz.col, imp.frac.min, imp.frac.max, 
                         veg.frac.min, veg.frac.max, tree.frac.min, tree.frac.max, 
                         stringsAsFactors = F)
    return(lcz.df)
}



# ---------------------------------------------------------------------------- #
# get FLUXNET tower info from file 'FLX_AA-Flx_BIF_LATEST.xlsx'
# ---------------------------------------------------------------------------- #
get.fluxnet.info <- function(file = '<path>/FLX_AA-Flx_BIF_LATEST.xlsx') {

    library(openxlsx)
    site.info <- openxlsx::read.xlsx(file) %>% 
                 filter(VARIABLE %in% c('COUNTRY', 'SITE_NAME', 'IGBP', 
                                        'LOCATION_LAT', 'LOCATION_LONG'))
    
    # three US sites have multiple locations, use the last rows of lat/lon 
    out.sites <- names(which(table(site.info$SITE_ID) > 5))
    out.info <- site.info %>% filter(SITE_ID %in% out.sites) %>% 
                              group_by(SITE_ID) %>%
                              filter(row_number() %in% c(1:3, n() - 1, n())) %>% 
                              ungroup()

    all.info <- rbind(site.info %>% filter(!SITE_ID %in% out.sites), out.info)

    ref.info <- all.info %>% reshape2::dcast(SITE_ID ~ VARIABLE, 
                                             value.var = 'DATAVALUE') %>% 
                             mutate(latitude = as.numeric(LOCATION_LAT), 
                                    longitude = as.numeric(LOCATION_LONG), 

                                    # change WSA, OSH, CSH to WSAV, OSHR, CSHR
                                    IGBP = ifelse(IGBP == 'WSA', 'WSAV', IGBP), 
                                    IGBP = ifelse(IGBP == 'OSH', 'OSHR', IGBP),
                                    IGBP = ifelse(IGBP == 'CSH', 'CSHR', IGBP))
    
    return(ref.info)
}   # end of get.fluxnet.info()




# ---------------------------------------------------------------------------- #
# function that store the biomes-specific RMSE of predicted vs observed Reco
# values here are consistent with Wu et al. (in prep) given presented NN model
# Dien Wu, last update on 01/03/2020
# RMSE is calculated from NN_train_reco_glb.r
# ---------------------------------------------------------------------------- #
get.reco.cv <- function(train.feature = c('fluxnet', 'era5', 'daymet_nldas')[1], 
                        reg = c('conus', 'global')) {
    
    biomes <- c('CRO', 'CRO/NVM',  'CSHR',      'DBF',    'DNF', 
                'EBF',     'ENF',   'GRA',       'MF',   'OSHR',   
                'SAV',     'WET',  'WSAV',  'ENF/DBF', 'EBF/DBF')

    # if using NN models trained by global flux sites
    if (reg == 'global') {
        if (train.feature == 'era5') { 
            cv <- c(0.457,    0.457,    0.272,    0.421,  0.871,  
                    0.399,    0.423,    0.534,    0.442,  0.715,   
                    0.361,    0.433,    0.353,    0.422,  0.406)
            # trained by 93 global sites (data 2010-2014)

        } else if (train.feature == 'fluxnet') {
            cv <- c(0.362,    0.380,    0.250,  0.344,   0.955, 
                    0.351,    0.386,    0.428,  0.429,   0.789, 
                    0.291,    0.307,    0.344,  0.372,   0.359)  
            # trained by 93 global sites (data 2010-2014)

        } else {
            stop(paste('get.reco.cv(): NO Reco uncertainty found for training feature from', 
                        train.feature))
        }   # end if 
    }   # end if reg GLOBAL

    if (reg == 'conus') {
        if (train.feature == 'era5') { 
            cv <- c(0.441,    0.441,       NA,    0.361,     NA,  
                       NA,    0.405,    0.496,    0.633,  0.895,   
                       NA,    0.353,    0.514,    0.384,  0.361)
            # trained by 27 US sites (data 2010-2014)

        } else if (train.feature == 'fluxnet') {
            cv <- c(0.404,    0.404,       NA,    0.315,      NA, 
                       NA,    0.368,    0.473,    0.612,   0.434, 
                       NA,    0.287,    0.498,    0.343,   0.315)  
            # trained by 27 US sites (data 2010-2014)

        } else if (train.feature == 'daymet_nldas') {
            cv <- c(0.425,    0.425,       NA,    0.340,     NA, 
                       NA,    0.411,    0.495,    0.667,  0.628,
                       NA,    0.344,    0.490,    0.377,  0.315)
            # trained by 27 US sites (data 2010-2014)
        }
    }   # end if reg CONUS
   

    cv.df <- data.frame(biomes, cv, stringsAsFactors = F)
    return(cv.df)
    
}   # end of get.reco.cov()



# ---------------------------------------------------------------------------- #
# function to find all 4-day or 8-day date strings given timestr (as YYYYMMDD), 
# and nhrs (hours backwards or forwards), 
# needed for get CSIF or GOSIF or MODIS data, DW, 04/09/2019

# nday: day interval, e.g., every 4, 8, 16 days
# timestr: in the form of YYYYMMDDHH in UTC, HH for overpass hour
# nhrs: how many hours backwards, if NULL, no need to calculate `start.*`
# ---------------------------------------------------------------------------- #
get.doy.from.timestr <- function(timestr, nday = c(4, 8, 16)[1], nhrs = NULL, 
                                 tz = 'UTC') {
    
    days <- seq(1, 366, nday)

    # get time info for input `timestr`
    yr0 <- as.numeric(substr(timestr, 1, 4))
    hr0 <- as.numeric(substr(timestr, 9, 10))
    date0  <- as.POSIXlt(as.character(timestr), format = '%Y%m%d%H', tz = tz)
    yrdoy0 <- as.numeric(strftime(date0, format = '%Y%j', tz = tz))
    doy0   <- as.numeric(substr(yrdoy0, 5, 7))

    # find the DOY and timestr for the nearest day interval
    find.doy   <- days[findInterval(doy0, days)]
    find.yrdoy <- as.numeric(paste0(yr0, formatC(find.doy, width = 3, flag = 0)))
    find.date  <- as.Date(find.doy - 1, origin = paste0(yr0, '-01-01'))
    find.timestr <- paste0(strftime(find.date, '%Y%m%d', tz = tz), 
                           formatC(hr0, width = 2, flag = 0))
    
    # create time string info
    time.info <- data.frame(timestr0 = timestr, date0 = date0, 
                            yr0 = yr0, doy0 = doy0, yrdoy0 = yrdoy0, 
                            find.timestr = find.timestr, find.date = find.date, 
                            find.doy = find.doy, find.yrdoy = find.yrdoy, 
                            stringsAsFactors = F)


    ## if nhrs is given, need to trace backwards or forwards by nhrs, 
    if (!is.null(nhrs)) {        # nhrs can be positive or negative
        
        date2 <- date0 + nhrs * 60 * 60
        timestr2 <- strftime(date2, '%Y%m%d%H', tz = tz)
        yr2    <- as.numeric(substr(timestr2, 1, 4))
        yrdoy2 <- as.numeric(strftime(date2, format = '%Y%j', tz = tz))
        doy2   <- as.numeric(substr(yrdoy2, 5, 7))
        find.doy2 <- days[findInterval(doy2, days)]

        time.info <- time.info %>% 
                     mutate(nhrs.timestr = timestr2, nhrs.date = date2, 
                            nhrs.yr = yr2, nhrs.doy = doy2, nhrs.yrdoy = yrdoy2, 
                            min.timestr = min(timestr0, timestr2, find.timestr), 
                            max.timestr = max(timestr0, timestr2, find.timestr), 
                            min.doy = min(doy0, find.doy, doy2), 
                            max.doy = max(doy0, find.doy, doy2))

                        
    } else {
        time.info <- time.info %>% 
                     mutate(min.timestr = min(timestr0, find.timestr), 
                            max.timestr = max(timestr0, find.timestr), 
                            min.doy = min(doy0, find.doy), 
                            max.doy = max(doy0, find.doy))

    } # end if 
    

    return(time.info)
}

   