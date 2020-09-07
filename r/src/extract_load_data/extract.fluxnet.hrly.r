#' subroutine to get daily mean fluxes from FLUXNET, DW
#' extract observed GPP/reco/NEE for every 4 days (according to modeled GPP resoulation)
#' this function requires grab.fluxnet.obs()

#' @param flux.xlsx.file filepath for a file called 'FLX_AA-Flx_BIF_LATEST.xlsx', 
#'                      when downloaded the data from FLUXNET2015 

#' @param site.name site name, e.g., 'US-NR1'
#' @param flux.path path that store raw csv for @param site.name
#'                   '<path>/FLX_**-***_FLUXNET2015_FULLSET_HH_YYYY-YYY_*-3.csv'

# --------------------------------------------------------------------------- #
extract.flux.hrly <- function(flux.xlsx.file, site.name = 'US-NR1', flux.path, 
                              outpath, overwriteTF = F, workdir) {

  try({
    
    # load dependencies 
    setwd(workdir); source(file.path(workdir, 'r/dependencies.r'), local = T)

    # get info from FLUXNET folder names
    flux.files <- list.files(flux.path, 'FULLSET_H', recursive = T, full.names = T)
    year.info <- strsplit.to.df(basename(flux.files), sep = '_') %>% 
                 mutate(min.year = as.numeric(substr(V6, 1, 4)), 
                        max.year = as.numeric(substr(V6, 6, 10)), 
                        site = V2, country = substr(V2, 1, 2), 
                        path = flux.files) %>% dplyr::select(-starts_with('V')) 

    # get site info, e.g., lat, lon, igbp and time zone
    all.info <- get.fluxnet.info(flux.xlsx.file) %>% 
                right_join(year.info, by = c('SITE_ID' = 'site')) %>%
                mutate(tz = tz_lookup_coords(latitude, longitude, method = 'accurate')) %>% 
                dplyr::rename(site = SITE_ID, country.name = COUNTRY, full.name = SITE_NAME)
    
    # ----------------------------------------------------------------------- #
    cat(paste('\n\n# ---- working on', site.name, ' ---- #\n'))
    sel.info <- all.info %>% filter(site == site.name)
    min.yr <- sel.info$min.year
    max.yr <- sel.info$max.year 
    igbp   <- sel.info$IGBP

    # ------------------------------------------------------------------- #
    obs.fn <- file.path(outpath, paste0('hrly_OBS_FLUXNET_', site.name, '_', igbp, 
                                        '_', min.yr, '_', max.yr, '.rds'))
    print(obs.fn)

    if ( file.exists(obs.fn) & overwriteTF == F ) {
        obs.list <- readRDS(obs.fn)

    } else {

        # call grab.fluxnet.obs(), see script below
        obs.list <- grab.fluxnet.hrly.obs(file = sel.info$path, min.yr, max.yr, 
                                          time.zone = sel.info$tz) 
        
        # add useful site info
        obs.list$mins.dat <- obs.list$mins.dat %>% 
                             mutate(site.name = site.name, 
                                    site.igbp = sel.info$IGBP, 
                                    site.lati = sel.info$latitude, 
                                    site.long = sel.info$longitude)

        obs.list$hr.dat <- obs.list$hr.dat %>% 
                           mutate(site.name = site.name, site.igbp = sel.info$IGBP, 
                                  site.lati = sel.info$latitude, 
                                  site.long = sel.info$longitude)

        if (!is.null(obs.list)) {
            obs.list <- c(sel.info, obs.list)
            saveRDS(obs.list, obs.fn) 
        }  # save data in rds 
    }   # end if
  
  })    # end of try()

} # end of subroutine 



# --------------------------------------------------------------------------- #
# inner subroutine to grab observed FLUXNET fluxes, DW, 05/09/2019 
# output fluxes: 30mins, hourly and daily data
# bug fixed when averaging fluxes, remove data whose QC == 3, DW, 05/31/2019
# no need to average to daily or month scales, DW, 08/09/2020
# time zone of hourly fluxes is local standard time, DW, 08/16/2020

# ** Time zone convention, https://fluxnet.org/data/aboutdata/data-variables/
# Time must the reported in local standard time (i.e., without “Daylight Saving Time”). 
# The time zone must be specified using the BADM template for the site.

grab.fluxnet.hrly.obs <- function(file, min.yr = NULL, max.yr = NULL, time.zone) {

    # ----------------------------------------------------------------------- #
    # read initial data, every 30mins fluxes
    cat('grab.fluxnet.hrly.obs(): reading inital 30mins GPP fluxes...it takes time...\n')
    dat <- read.csv(file, sep = ',', header = T)

    # select gap-filled columns
    init.dat <- dat %>% dplyr::select(starts_with('TIMESTAMP'), ends_with('MEAN'), 
                                      starts_with('TA_F'), starts_with('TS_F'),
                                      starts_with('SWC_F'), starts_with('VPD_F'),
                                      starts_with('RH'), ends_with('_MEAN_QC'), 
                                      ends_with('_REF_QC'), ends_with('_REF'))  
    
    print(colnames(init.dat))
    cat('REMOVING BAD DATA (QC == 3), QC range before and after: \n')
    print(range(init.dat$NEE_VUT_REF_QC))

    # remove all poor gapfilled data, QC = 3 and missing data of -9999
    init.dat <- init.dat %>% mutate(timestr.lt = as.character(TIMESTAMP_START)) %>% 
                             filter_at(vars(ends_with('QC')), all_vars(.!= 3)) %>% 
                             filter_at(vars(ends_with('MEAN')), all_vars(.!= -9999)) 

    # select data within min.yr and max.yr, 
    if (!is.null(min.yr)) init.dat <- init.dat %>% filter(timestr.lt >= paste0(min.yr, '01010000')) 
    if (!is.null(max.yr)) init.dat <- init.dat %>% filter(timestr.lt <= paste0(max.yr, '12310000'))
    if (nrow(init.dat) == 0) stop('NO obs data with fine quality (QC != 3) between the desired time period')
    print(range(init.dat$NEE_VUT_REF_QC))
    #debug <- init.dat %>% filter(NEE_VUT_REF_QC == 0)


    # debug daylight saving conversion, force time zone as standard time 
    # since FLUXNET always use standard time without daylight saving, DW, 08/22/2020
    # by default, R will use both MDT and MST if applied
    tmp.dat <- init.dat %>% filter(substr(TIMESTAMP_START, 1, 4) == '2010')
    tmp.dat$date.lt = as.POSIXct(tmp.dat$timestr.lt, tz = time.zone, format = '%Y%m%d%H%M')
    tmp.dat$date.lt[1E4:10010]
    tmp.dat$date.lt[1:10]

    init.dat$date.utc = init.dat$date.lt
    attributes(init.dat$date.utc)$tzone <- 'UTC'
    init.dat$timestr.utc <- format(init.dat$date.utc, format = '%Y%m%d%H%M')

    # add time string, and convert local STANDARD time to UTC
    # adjust daytime saving if there is
    init.dat$date.lt = as.POSIXct(init.dat$timestr.lt, tz = time.zone, format = '%Y%m%d%H%M')
    init.dat$date.utc = init.dat$date.lt
    attributes(init.dat$date.utc)$tzone <- 'UTC'
    init.dat$timestr.utc <- format(init.dat$date.utc, format = '%Y%m%d%H%M')
    #debug <- init.dat %>% group_by(substr(timestr.utc, 9, 10)) %>% summarise_all(mean) 


    # ----------------------------------------------------------------------- #
    # average all variables from every 30 mins to every hour
    cat('grab.fluxnet.hrly.obs(): calculating hourly fluxes and variables; it takes time...\n')

    # replace all -9999 with NA, will be removed afterwards
    #hr.dat <- init.dat %>% dplyr::select(-c('TIMESTAMP_START', 'TIMESTAMP_END'))
    hr.dat <- init.dat %>% mutate(timestr.utc = substr(timestr.utc, 1, 10), 
                                  timestr.lt = substr(timestr.lt, 1, 10)) 
    hr.dat[hr.dat == -9999] <- NA   

    # abandon the 30mins interval where only one time slot is available, DW, 08/16/2020
    if (!grepl('FULLSET_HR_', file)) {
        rm.mins <- names( which(table(hr.dat$timestr.utc) == 1))
        hr.dat <- hr.dat %>% filter(!timestr.utc %in% rm.mins) 
    }

    hr.dat <- hr.dat %>% group_by(timestr.utc, timestr.lt) %>% 
              summarise_if(is.numeric, mean, na.rm = T) %>% ungroup() %>% 
              mutate(hr.utc = as.numeric(substr(timestr.utc, 9, 10)), 
                     hr.lt = as.numeric(substr(timestr.lt, 9, 10)))
    
    #debug <- hr.dat %>% group_by(hr.utc) %>% summarise_all(mean) 

    obs.list <- list(mins.dat = init.dat, hr.dat = hr.dat)
    return(obs.list)
}

# end of subroutine





# --------------------------------------------------------------------------- #
# inner subroutine to grab observed AmeriFlux fluxes, DW, 10/14/2019 
# output fluxes: 30mins, hourly and daily data
grab.ameriflux.hrly.obs <- function(file, site.name, site.igbp, time.zone = 'UTC', 
                                    min.yr = NULL, max.yr = NULL) {

    library(lubridate)

    # ----------------------------------------------------------------------- #
    # read initial data, every 30mins fluxes
    cat('grab.ameriflux.hrly.obs(): reading inital 30mins GPP fluxes...it takes time...\n')
    dat <- read.csv(file, sep = ',', header = T, skip = 2)
    
    # select gap-filled columns
    init.dat <- dat %>% dplyr::select(starts_with('TIMESTAMP'), ends_with('CO2'), 
                                      starts_with('TA'), starts_with('TS'),
                                      starts_with('SWC'), starts_with('VPD'), 
                                      starts_with('GPP'), starts_with('RECO')) 
    print(colnames(init.dat))

    # select data within min.yr and max.yr, 
    # remove all poor gapfilled data, QC = 3 and missing data of -9999
    init.dat <- init.dat %>% mutate(timestr = as.character(TIMESTAMP_START)) #%>% 
                             #filter_at(vars(ends_with('QC')), all_vars(.!= 3))

    # select data within min.yr and max.yr, 
    if (!is.null(min.yr)) init.dat <- init.dat %>% filter(timestr >= paste0(min.yr, '01010000')) 
    if (!is.null(max.yr)) init.dat <- init.dat %>% filter(timestr <= paste0(max.yr, '01010000'))
    if (nrow(init.dat) == 0) stop('NO obs data between the desired time period')
    init.dat[init.dat == -9999] <- NA

    # add time string
    init.dat <- init.dat %>% mutate(site.name = site.name, igbp = site.igbp,
                                    hr = as.numeric(substr(timestr, 9, 10)), 
                                    LT = as.POSIXct(timestr, format = '%Y%m%d%H', tz = time.zone), 
                                    DoY = lubridate::yday(LT)) #%>% 
                                    
    # ----------------------------------------------------------------------- #
    # average all variables from every 30 mins to every hour
    cat('grab.ameriflux.hrly.obs(): calculating hourly fluxes and variables; it takes time...\n')

    # replace all -9999 with NA, will be removed afterwards
    hr.dat <- init.dat %>% dplyr::select(-c('TIMESTAMP_START', 'TIMESTAMP_END'))
    hr.dat[hr.dat == -9999] <- NA   
    hr.dat <- hr.dat %>% mutate(timestr = substr(timestr, 1, 10)) %>% 
              group_by(timestr) %>% summarise_if(is.numeric, mean, na.rm = T) %>% 
              ungroup() %>% mutate(LT = as.POSIXct(timestr, format = '%Y%m%d%H', tz = time.zone), 
                                   DoY = yday(LT), hr = as.numeric(substr(timestr, 9, 10)), 
                                   site.name = site.name, igbp = site.igbp)
              
    # average GPP from hourly to daily
    cat('grab.ameriflux.hrly.obs(): calculating daily fluxes and variables; it takes time...\n')
    day.dat <- hr.dat %>% dplyr::select(-'hr') %>% 
               mutate(timestr = substr(timestr, 1, 8)) %>%
               group_by(timestr) %>% summarize_if(is.numeric, mean, na.rm = T) %>% 
               ungroup() %>% mutate(LT = as.POSIXct(timestr, format = '%Y%m%d', tz = time.zone), 
                                    DoY = yday(LT), site.name = site.name, 
                                    igbp = site.igbp)
               
    obs.list <- list(mins.dat = init.dat, hr.dat = hr.dat, day.dat = day.dat)
    return(obs.list)
}
