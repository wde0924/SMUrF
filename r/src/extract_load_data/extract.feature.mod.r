#' function to extract tair, tsoil, CSIF, and GPP from models at flux sites
#' written by Dien Wu, 07/30/2019 

#' @param site.name unique site name

extract.feature.mod <- function(outpath, 
                                site.name = 'US-NR1', 
                                flux.source = c('FLUXNET', 'AMERIFLUX')[1], 

                                # ERA5, NLDAS and Daymet paths
                                # and their available years
                                nldas.path = NULL, daymet.path = NULL, era5.path, 
                                #.min.yr = 2010, t.max.yr = 2018, 
                                
                                # SIF info
                                csif.path, csif.min.yr = 2010, 
                                tsif.path = NULL, 

                                # path for SMUrF repo and tmp directory
                                smurf_wd, tmpdir) {

    # Ensure dependencies are loaded for current node/process
    setwd(smurf_wd); source('r/dependencies.r', local = T)
    if (!is.null(tmpdir)) raster::rasterOptions(tmpdir = tmpdir)
    nhrs <- 23 

    # ----------------------------------------------------------------------- #
    if (flux.source == 'FLUXNET') {
        
        rds.files <- list.files(outpath, 'daily_OBS_FLUXNET_', full.names = T)
        rds.files <- rds.files[!grepl('model', rds.files)]
        site.info.all <- strsplit.to.df(gsub('.rds', '',basename(rds.files))) %>%
                         mutate(min.year = as.numeric(V6), max.year = as.numeric(V7),
                                site = V4, igbp = V5, country = substr(site, 1, 2), 
                                rds.file = rds.files) %>% dplyr::select(-starts_with('V')) 

        # read FLUXNET data points from .rds file
        sel.info  <- site.info.all %>% filter(site == site.name)
        daily.obs <- readRDS(sel.info$rds.file)$daily.obs

    } else if (flux.source == 'AMERIFLUX') {
        
        rds.files <- list.files(outpath, 'hrly_OBS_AMERIFLUX_', full.names = T)
        rds.files <- rds.files[!grepl('model', rds.files)]
        site.info.all <- strsplit.to.df(gsub('.rds', '',basename(rds.files))) %>%
                         mutate(min.year = as.numeric(V6), max.year = as.numeric(V7),
                                site = V4, igbp = V5, country = substr(site, 1, 2), 
                                rds.file = rds.files) %>% dplyr::select(-starts_with('V')) 
        
        # read AMERIFLUX data points from .rds file
        sel.info  <- site.info.all %>% filter(site == site.name)
        daily.obs <- readRDS(sel.info$rds.file)$day.dat %>% na.omit() 
    }   # end if

    sel.obs <- daily.obs %>% filter(substr(timestr, 1, 4) >= csif.min.yr) %>% 
               dplyr::rename(site = site.name)
    if (nrow(sel.obs) == 0) {cat('NO data beyond CSIF years...\n'); return()}

    # **** site.nameal extent, deterimined by sites
    site.lat <- readRDS(sel.info$rds.file)$site.info$latitude 
    site.lon <- readRDS(sel.info$rds.file)$site.info$longitude 
    site.loc <- data.frame(longitude = site.lon, latitude = site.lat)

    ext <- extent(floor(site.lon - 0.5), ceiling(site.lon + 0.5), 
                  floor(site.lat - 0.5),  ceiling(site.lat + 0.5))

    # ----------------------------------------------------------------------- #
    if (sel.info$min.year < csif.min.yr) sel.info$min.year <- csif.min.yr 
    if (flux.source == 'FLUXNET') head.string <- 'daily_OBS_FLUXNET_model_'
    if (flux.source == 'AMERIFLUX') head.string <- 'daily_OBS_AMERIFLUX_model_'

    rds.fn <- file.path(outpath, paste0(head.string, sel.info$site, '_', 
                                        sel.info$igbp, '_', sel.info$min.year, 
                                        '_', sel.info$max.year, '.rds')) 
    print(rds.fn)

    ## get unique time stamps, site names, and lon/lat 
    uni.timestr <- unique(sel.obs$timestr)  # YYYYMMDD
    if ('20180307' %in% uni.timestr) sel.obs$TSIF_inst <- NA 
    
    for (t in 1 : length(uni.timestr)) {  # Jan 1, Mar 59, June 151, Sept 241
        
        tmp.timestr <- paste0(uni.timestr[t], '00')
        tmp.yr <- substr(tmp.timestr, 1, 4)
        if (t %% 10 == 0) 
            cat(paste('\n# --- Working on site', site.name, '; time:', tmp.timestr,
                       ';', signif(t / length(uni.timestr)) * 100, '% done --- #\n'))

        #if (tmp.yr >= t.min.yr & tmp.yr <= t.max.yr) {
            
        # calculate daily mean Tair and Tsoil, nhrs = 23 -> 00:00 to 23:00, one day
        # e.g., nhrs = 23; 23 hours forward for calculate daily mean temp
        # soil temp and water content, both from Level 1, the topest layer
        ta.era5.brk  <- prep.era5(era5.path, '2T', tmp.timestr, ext, nhrs)
        ts.era5.brk  <- prep.era5(era5.path, 'STL1', tmp.timestr, ext, nhrs)
        swv.era5.brk <- prep.era5(era5.path, 'SWVL1', tmp.timestr, ext, nhrs)
        
        ta.era5.rt  <- mean(ta.era5.brk)    # degC
        ts.era5.rt  <- mean(ts.era5.brk)    # degC
        swv.era5.rt <- mean(swv.era5.brk)   # m3/m3, unitless

        ## grab the closest 4-day CSIF and daily TROPOMI SIF, they can be negative
        csif.mean.rt <- grab.csif(csif.path, tmp.timestr, ext, var = 'clear_daily_SIF') 
        csif.inst.rt <- grab.csif(csif.path, tmp.timestr, ext, var = 'clear_inst_SIF') 

        # extract CSIF and ERA5 based temps, default method is 'simple'
        sel.obs[t, 'CSIF_mean'] <- raster::extract(csif.mean.rt, site.loc)
        sel.obs[t, 'CSIF_inst'] <- raster::extract(csif.inst.rt, site.loc)
        sel.obs[t, 'TA_era5']  <- raster::extract(ta.era5.rt, site.loc)
        sel.obs[t, 'TS_era5']  <- raster::extract(ts.era5.rt, site.loc)
        sel.obs[t, 'SWV_era5'] <- raster::extract(swv.era5.rt, site.loc)

        # additional extractions for US sites from daymet and NLDAS
        if (sel.info$country == 'US') {
            
            ts.nldas.brk  <- prep.soil.nldas(nldas.path, 'TSOIL', tmp.timestr, ext, nhrs)
            swc.nldas.brk <- prep.soil.nldas(nldas.path, 'SOILM', tmp.timestr, ext, nhrs)
            ta.daymet.brk <- prep.tair.daymet(daymet.path, timestr = tmp.timestr, ext = ext, nhrs = nhrs)

            ts.nldas.rt  <- mean(ts.nldas.brk)   # degC
            swc.nldas.rt <- mean(swc.nldas.brk)  # unit in kg/m2
            ta.daymet.rt <- mean(ta.daymet.brk)  # degC

            # if during TROPOMI years, grab 500m TROPOMI SIF 
            tsif.inst.rt <- grab.tsif(tsif.path, tmp.timestr, ext, var = 'sif')

            if (!is.null(tsif.inst.rt)) 
                sel.obs[t, 'TSIF_inst'] <- raster::extract(tsif.inst.rt, site.loc)
                
            sel.obs[t, 'TA_daymet'] <- raster::extract(ta.daymet.rt, site.loc)
            sel.obs[t, 'TS_nldas']  <- raster::extract(ts.nldas.rt, site.loc)
            sel.obs[t, 'SWC_nldas'] <- raster::extract(swc.nldas.rt, site.loc)
        }  # end if US
    
        #}   # end if tmp.yr
      
    }   # end of t

    saveRDS(sel.obs, file = rds.fn)

}   # end of subroutinr



# --------------------------------------------------------------------------- #
# interpolate SIF for CONUS
extract.feature.sif.amer <- function(outpath, site.name = 'US-NR1',
                                     csif.path, tsif.path, smurf_wd, tmpdir) {

    # Ensure dependencies are loaded for current node/process
    setwd(smurf_wd); source('r/dependencies.r', local = T)
    if (!is.null(tmpdir)) raster::rasterOptions(tmpdir = tmpdir)

    rds.files <- list.files(outpath, 'hrly_OBS_AMERIFLUX_', full.names = T)
    rds.files <- rds.files[!grepl('model', rds.files)]
    site.info.all <- strsplit.to.df(gsub('.rds', '',basename(rds.files))) %>%
                     mutate(min.year = as.numeric(V6), max.year = as.numeric(V7),
                            site = V4, igbp = V5, country = substr(site, 1, 2), 
                            rds.file = rds.files) %>% dplyr::select(-starts_with('V')) 
    
    # read AMERIFLUX data points from .rds file
    sel.info  <- site.info.all %>% filter(site == site.name)
    hr.obs <- readRDS(sel.info$rds.file)$hr.dat

    # only select observations at 1pm local time
    sel.obs <- hr.obs %>% filter(timestr >= '2018060100', timestr < '2018090100',
                                 !is.na(GPP_PI_F), hr == 13) %>% 
                          dplyr::rename(site = site.name)
    print(paste(site.name, nrow(sel.obs), 'rows'))

    # **** site.nameal extent, deterimined by sites
    site.lat <- readRDS(sel.info$rds.file)$site.info$latitude 
    site.lon <- readRDS(sel.info$rds.file)$site.info$longitude 
    site.loc <- data.frame(longitude = site.lon, latitude = site.lat)
    ext <- extent(floor(site.lon - 0.1), ceiling(site.lon + 0.1), 
                  floor(site.lat - 0.1),  ceiling(site.lat + 0.1))

    # ----------------------------------------------------------------------- #
    rds.fn <- file.path(outpath, paste0('hrly_OBS_AMERIFLUX_SIF_', sel.info$site, 
                                        '_', sel.info$igbp, '_JJA2018.rds')) 
    print(rds.fn)

    ## get unique time stamps, site names, and lon/lat 
    uni.timestr <- unique(sel.obs$timestr)  # YYYYMMDD
    for (t in 1 : length(uni.timestr)) {  # Jan 1, Mar 59, June 151, Sept 241
        
        tmp.timestr <- uni.timestr[t]
        tmp.yr <- as.numeric(substr(tmp.timestr, 1, 4))

        cat(paste('\n\n# --- Working on site', site.name, '; time:', tmp.timestr,
                    ';', signif(t / length(uni.timestr)) * 100, '% done --- #\n'))
    
        # if during TROPOMI years, grab 500m TROPOMI SIF 
        source(file.path(smurf_wd, 'r/dependencies.r'))
        tsif.inst.df <- grab.tsif(tsif.path, tmp.timestr, minlon = site.lon - 0.05, 
                                  maxlon = site.lon + 0.05, minlat = site.lat - 0.05, 
                                  maxlat = site.lat + 0.05)
        #tsif.inst.rt <- df2raster(tsif.inst.df[,c('lon', 'lat', 'SIF')])

        uni.lat <- unique(sort(tsif.inst.df$lat))
        uni.lon <- unique(sort(tsif.inst.df$lon))
        lat.indx <- findInterval(site.lat, uni.lat)
        lon.indx <- findInterval(site.lon, uni.lon)
        tsif.site <- tsif.inst.df %>% filter(lon == uni.lon[lon.indx], lat == uni.lat[lat.indx])
        sel.obs[t, 'TSIF_inst'] <- tsif.site$SIF

        csif.inst.rt <- grab.csif(csif.path, tmp.timestr, ext, var = 'clear_inst_SIF') 
        sel.obs[t, 'CSIF_inst'] <- raster::extract(csif.inst.rt, site.loc)
        print(sel.obs[t, c('GPP_PI_F', 'TSIF_inst', 'CSIF_inst')]) 

    }   # end of t

    saveRDS(sel.obs, file = rds.fn)
}   # end of function 


# --------------------------------------------------------------------------- #
# interpolate temps
extract.feature.mod.amer <- function(outpath, site.name = 'US-NR1',
                                     nldas.path, daymet.path, era5.path, 
                                     csif.path, tsif.path, smurf_wd, tmpdir) {

    # Ensure dependencies are loaded for current node/process
    setwd(smurf_wd); source('r/dependencies.r', local = T)
    if (!is.null(tmpdir)) raster::rasterOptions(tmpdir = tmpdir)

    rds.files <- list.files(outpath, 'hrly_OBS_AMERIFLUX_', full.names = T)
    rds.files <- rds.files[!grepl('model', rds.files)]
    site.info.all <- strsplit.to.df(gsub('.rds', '',basename(rds.files))) %>%
                     mutate(min.year = as.numeric(V6), max.year = as.numeric(V7),
                            site = V4, igbp = V5, country = substr(site, 1, 2), 
                            rds.file = rds.files) %>% dplyr::select(-starts_with('V')) 
    
    # read AMERIFLUX data points from .rds file
    sel.info  <- site.info.all %>% filter(site == site.name)
    hr.obs <- readRDS(sel.info$rds.file)$hr.dat

    # only select observations at 1pm local time
    sel.obs <- hr.obs %>% filter(timestr >= '2018060100', timestr < '2018090100',
                                 !is.na(GPP_PI_F), hr == 13) %>% 
                          dplyr::rename(site = site.name)
    print(paste(site.name, nrow(sel.obs), 'rows'))

    # **** site.nameal extent, deterimined by sites
    site.lat <- readRDS(sel.info$rds.file)$site.info$latitude 
    site.lon <- readRDS(sel.info$rds.file)$site.info$longitude 
    site.loc <- data.frame(longitude = site.lon, latitude = site.lat)

    ext <- extent(floor(site.lon - 0.5), ceiling(site.lon + 0.5), 
                  floor(site.lat - 0.5),  ceiling(site.lat + 0.5))

    # ----------------------------------------------------------------------- #
    rds.fn <- file.path(outpath, paste0('daily_OBS_AMERIFLUX_model_', sel.info$site, 
                                        '_', sel.info$igbp, '_JJA2018.rds')) 
    print(rds.fn)

    ## get unique time stamps, site names, and lon/lat 
    uni.timestr <- unique(sel.obs$timestr)  # YYYYMMDD

    # e.g., nhrs = 23; 23 hours forward for calculate daily mean temp
    # e.g., nhrs = NULL, instant temps
    nhrs <- NULL 

    for (t in 1 : length(uni.timestr)) {  # Jan 1, Mar 59, June 151, Sept 241
        
        tmp.timestr <- uni.timestr[t]
        tmp.yr <- as.numeric(substr(tmp.timestr, 1, 4))

        cat(paste('\n\n# --- Working on site', site.name, '; time:', tmp.timestr,
                    ';', signif(t / length(uni.timestr)) * 100, '% done --- #\n'))
    
        # if during TROPOMI years, grab 500m TROPOMI SIF 
        source(file.path(smurf_wd, 'r/dependencies.r'))
        tsif.inst.df <- grab.tsif(tsif.path, tmp.timestr, minlon = site.lon - 0.05, 
                                  maxlon = site.lon + 0.05, minlat = site.lat - 0.05, 
                                  maxlat = site.lat + 0.05)
        #tsif.inst.rt <- df2raster(tsif.inst.df[,c('lon', 'lat', 'SIF')])

        uni.lat <- unique(sort(tsif.inst.df$lat))
        uni.lon <- unique(sort(tsif.inst.df$lon))
        lat.indx <- findInterval(site.lat, uni.lat)
        lon.indx <- findInterval(site.lon, uni.lon)
        tsif.site <- tsif.inst.df %>% filter(lon == uni.lon[lon.indx], lat == uni.lat[lat.indx])
        sel.obs[t, 'TSIF_inst'] <- tsif.site$SIF

        # calculate daily mean Tair and Tsoil, nhrs = 23 -> 00:00 to 23:00, one day
        # soil temp and water content, both from Level 1, the topest layer
        ta.era5.brk  <- prep.era5(era5.path, '2T', tmp.timestr, ext, nhrs)
        ts.era5.brk  <- prep.era5(era5.path, 'STL1', tmp.timestr, ext, nhrs)
        swv.era5.brk <- prep.era5(era5.path, 'SWVL1', tmp.timestr, ext, nhrs)
        
        ta.era5.rt  <- mean(ta.era5.brk)    # degC
        ts.era5.rt  <- mean(ts.era5.brk)    # degC
        swv.era5.rt <- mean(swv.era5.brk)   # m3/m3, unitless

        ## grab the closest 4-day CSIF and daily TROPOMI SIF, they can be negative
        csif.mean.rt <- grab.csif(csif.path, tmp.timestr, ext, var = 'clear_daily_SIF') 
        csif.inst.rt <- grab.csif(csif.path, tmp.timestr, ext, var = 'clear_inst_SIF') 

        # extract CSIF and ERA5 based temps, default method is 'simple'
        sel.obs[t, 'CSIF_mean'] <- raster::extract(csif.mean.rt, site.loc)
        sel.obs[t, 'CSIF_inst'] <- raster::extract(csif.inst.rt, site.loc)
        sel.obs[t, 'TA_era5']  <- raster::extract(ta.era5.rt, site.loc)
        sel.obs[t, 'TS_era5']  <- raster::extract(ts.era5.rt, site.loc)
        sel.obs[t, 'SWV_era5'] <- raster::extract(swv.era5.rt, site.loc)

        # additional extractions for US sites from daymet and NLDAS
        ts.nldas.brk  <- prep.soil.nldas(nldas.path, 'TSOIL', tmp.timestr, ext, nhrs)
        swc.nldas.brk <- prep.soil.nldas(nldas.path, 'SOILM', tmp.timestr, ext, nhrs)
        ta.daymet.brk <- prep.tair.daymet(daymet.path, timestr = tmp.timestr, ext = ext, nhrs = nhrs)

        ts.nldas.rt  <- mean(ts.nldas.brk)   # degC
        swc.nldas.rt <- mean(swc.nldas.brk)  # unit in kg/m2
        ta.daymet.rt <- mean(ta.daymet.brk)  # degC

        sel.obs[t, 'TA_daymet'] <- raster::extract(ta.daymet.rt, site.loc)
        sel.obs[t, 'TS_nldas']  <- raster::extract(ts.nldas.rt, site.loc)
        sel.obs[t, 'SWC_nldas'] <- raster::extract(swc.nldas.rt, site.loc)

    }   # end of t

    saveRDS(sel.obs, file = rds.fn)

}   # end of subroutinr