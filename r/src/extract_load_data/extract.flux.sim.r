# subroutine to extract sim or obs GPP values at points from GPP fields
# DW, 06/12/2019 

# pts: two columns data frame or matrix, first for lon, second for lat
# extract: whether just cropped the intial rasterStack or extract values from pts
extract.sim.pts <- function(sim.files, sim.varname = c('SIF', 'GPP', 'Reco', 'NEE')[1], 
                            pts, dlat = 1, dlon = 1, rds.file) {
    
    library(raster)
    sim.ext <- extent(stack(sim.files[1], varname = paste0(sim.varname, '_mean')))

    # select towers that fall within GPP fields, only 9 sites over western US
    pts <- pts %>% filter(latitude >= sim.ext[3], latitude <= sim.ext[4],
                          longitude >= sim.ext[1], longitude <= sim.ext[2])

    print(pts$id)
    nsite <- nrow(pts)

    # break entire GPP fields into 1x1deg tiles to speed up the extrapolation process
    sim.lon <- trunc(seq(sim.ext[1], sim.ext[2], dlon))    
    sim.lat <- trunc(seq(sim.ext[3], sim.ext[4], dlat))
    group.pts <- pts %>% mutate(lon = sim.lon[findInterval(longitude, sim.lon)], 
                                lat = sim.lat[findInterval(latitude, sim.lat)], 
                                loc = paste(lon, lat))  # lower left corner

    uni.loc <- unique(group.pts$loc)    # unique lat/long locations
    uni.loc <- strsplit.to.df(uni.loc, sep = ' ', colnms = c('lon', 'lat')) %>% 
               mutate_all(funs(as.numeric), colnames(uni.loc))
    ngroup <- nrow(uni.loc) 
    cat(paste(ngroup, 'group(s)\n'))

    # ----------------------------------------------------------------------- #
    # first crop modeled fields and then extract values at points
    sim.df <- NULL
    for (g in 1 : ngroup) {

        cat(paste('# ----- extract.sim.pts(): Loading and cropping', sim.varname,
                  'for group', g, '----- #\n'))

        tmp.loc <- uni.loc[g, ]
        tmp.pts <- group.pts %>% filter(lon == tmp.loc$lon, lat == tmp.loc$lat) 
        tmp.ext <- extent(tmp.loc$lon, tmp.loc$lon + dlon, 
                          tmp.loc$lat, tmp.loc$lat + dlat)  # tmp cropping extent

        # read variable and its uncertainty, this takes few mins
        cat(paste('working on', basename(sim.files[1]), '\n'))
        sim.stk <- crop(stack(sim.files[1], varname = paste0(sim.varname, '_mean')), tmp.ext)
        #sim.sd.stk <- crop(stack(sim.files[1], varname = paste0(sim.varname, '_sd')), tmp.ext)

        if (length(sim.files) > 1) {
            for (r in 2 : length(sim.files)) {
                cat(paste('working on', basename(sim.files[r]), '\n'))
                sim.stk <- stack(sim.stk, crop(stack(sim.files[r], 
                                 varname = paste0(sim.varname, '_mean')), tmp.ext))
                #sim.sd.stk <- stack(sim.sd.stk, crop(stack(sim.files[r], varname = paste0(sim.varname, '_mean')), tmp.ext))
            }   # end for r
        }

        # -------------------------------------------------------------------- #
        # then, extrapolate the values for a latlong coordinate 
        # using two methods, simple vs. bilinear
        cat('extract.sim.pts(): extracting values at 0.05 grid cell from cropped 1x1deg field...\n')
        extract.pts <- tmp.pts %>% dplyr::select('longitude', 'latitude')

        sim <- raster::extract(sim.stk, extract.pts, 'simple') 
        #sim.sd <- raster::extract(sim.sd.stk, extract.pts, 'simple') 
        sim.int <- raster::extract(sim.stk, extract.pts, 'bilinear') 
        #sim.sd.int <- raster::extract(sim.sd.stk, extract.pts, 'bilinear') 

        sec <- gsub('X', '', colnames(sim))
        vec <- c(sim.varname, #paste0(sim.varname, '.sd'), 
                 paste0(sim.varname, '.int'))#, paste0(sim.varname, '.sd.int'))

        sim.array <- array(0, dim = c(dim(sim), 2), 
                            dimnames = list(1 : nrow(sim), sec, vec))

        sim.array[,,1] <- sim#; sim.array[,,2] <- sim.sd
        sim.array[,,2] <- sim.int#; sim.array[,,4] <- sim.sd.int 
        extract.pts <- extract.pts %>% mutate(rownum = 1:nrow(extract.pts), id = tmp.pts$id)

        # form array into data frame
        cat('extract.sim.pts(): forming into dataframe...\n')
        tmp.df <- as.data.frame.table(sim.array, stringsAsFactors = F) %>% 
                  rename(rownum = Var1, sec = Var2, fac = Var3, sim = Freq) %>% 
                  mutate(rownum = as.numeric(rownum), sec = as.numeric(sec)) %>%
                  left_join(extract.pts, by = 'rownum') %>% 
                  mutate(date = as.POSIXct(sec, origin = '1970-01-01 00:00:00', tz = 'UTC'), 
                         timestr = format(date, '%Y%m%d')) %>% dplyr::select(-rownum) 
        
        sim.df <- rbind(sim.df, tmp.df)
        gc()

    }   # end for g

    sim.df$fac <- as.factor(sim.df$fac)
    sim.df <- reshape2::dcast(sim.df, ... ~ fac, value.var = 'sim')
    saveRDS(sim.df, file = rds.file)

    return(sim.df)
}


# ---------------------------------------------------------------------------- #
# extract an area from rasterStack
extract.sim.area <- function(sim.files, extract.ext, 
                             varname = c('SIF_mean', 'SIF_sd', 'GPP_mean', 'GPP_sd')[3]) {
    library(raster)

    # read GPP and its uncertainty, this takes few mins
    cat(paste('extract.sim.area(): working on', basename(sim.files[1]), '\n'))
    sim.stk <- crop(stack(sim.files[1], varname = varname), extract.ext)

    if (length(sim.files) > 1) {
        for (r in 2 : length(sim.files)) {
            cat(paste('extract.sim.area(): working on', basename(sim.files[r]), '\n'))
            sim.stk <- stack(sim.stk, crop(stack(sim.files[r], varname = varname), 
                                           extract.ext))
        }   # end for r
    }   # end if

    return(sim.stk)
}  # end of subroutine


# ---------------------------------------------------------------------------- #
# extract Reco for an area
extract.sim.reco.area <- function(reco.files, extract.ext) {

    reco.stk <- crop(stack(reco.files), extract.ext)
    timestr  <- strsplit.to.df(names(reco.stk))$V4
    names(reco.stk) <- as.POSIXct(timestr, format = '%Y%m%d', tz = 'UTC')
    return(reco.stk)
}

# ---------------------------------------------------------------------------- #
# extract an area from rasterStack
extract.sim.sif.area <- function(sif.files, extract.ext) {
    
    library(raster)

    # read GPP and its uncertainty, this takes few mins
    cat(paste('working on', basename(sif.files[1]), '\n'))
    sif.stk <- crop(stack(sif.files[1]), extract.ext)

    for (r in 2 : length(sif.files)) {
        cat(paste('working on', basename(sif.files[r]), '\n'))
        sif.stk <- stack(sif.stk, crop(stack(sif.files[r]), extract.ext))
    }   # end for r

    return(sif.stk)

}  # end of subroutine

