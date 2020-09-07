#' subroutine to get daily mean fluxes from FLUXNET, DW
#' extract observed GPP/reco/NEE for every 4 days (according to modeled GPP resoulation)
#' this function requires grab.fluxnet.obs()

#' @param flux.xlsx.file filepath for a file called 'FLX_AA-Flx_BIF_LATEST.xlsx', 
#'                      when downloaded the data from FLUXNET2015 

#' @param site.name site name, e.g., 'US-NR1'
#' @param flux.path path that store raw csv for @param site.name
#'                   '<path>/FLX_**-***_FLUXNET2015_FULLSET_HH_YYYY-YYY_*-3.csv'

# --------------------------------------------------------------------------- #
extract.fluxnet.daily <- function(flux.xlsx.file, site.name = 'US-NR1', flux.path, 
                                  outpath, overwriteTF = F, smurf_wd) {

  try({
    
    # load dependencies 
    source(file.path(smurf_wd, 'r/dependencies.r'), local = T)
    
    # get info from FLUXNET folder names
    flux.files <- list.files(flux.path, 'FULLSET_D', recursive = T, full.names = T)
    year.info <- strsplit.to.df(basename(flux.files), sep = '_') %>% 
                 mutate(min.year = as.numeric(substr(V6, 1, 4)), 
                        max.year = as.numeric(substr(V6, 6, 10)), 
                        site = V2, country = substr(V2, 1, 2), path = flux.files) %>% 
                 dplyr::select(-starts_with('V')) 

    # get site info, e.g., lat, lon, igbp and time zone
    all.info <- get.fluxnet.info(flux.xlsx.file) %>% 
                right_join(year.info, by = c('SITE_ID' = 'site')) %>%
                mutate(tz = suppressWarnings(tz_lookup_coords(latitude, longitude))) %>% 
                dplyr::rename(site = SITE_ID, country.name = COUNTRY, full.name = SITE_NAME)
    
    # ----------------------------------------------------------------------- #
    cat(paste('\n\n# ---- working on', site.name, ' ---- #\n'))
    sel.info <- all.info %>% filter(site == site.name)
    min.yr <- sel.info$min.year
    max.yr <- sel.info$max.year 
    igbp   <- sel.info$IGBP

    # ------------------------------------------------------------------- #
    obs.fn <- file.path(outpath, paste0('daily_OBS_FLUXNET_', site.name, '_', 
                                        igbp, '_', min.yr, '_', max.yr, '.rds'))
    print(obs.fn)

    if ( file.exists(obs.fn) & overwriteTF == F ) {
        daily.obs <- readRDS(obs.fn)

    } else {

        # call grab.fluxnet.obs(), see script below
        daily.obs <- grab.fluxnet.daily.obs(file = sel.info$path, site.name, 
                                            time.zone = sel.info$tz, min.yr, max.yr, 
                                            igbp = sel.info$IGBP)

        if (!is.null(daily.obs)) {
            daily.obs <- list(site.info = sel.info, daily.obs = daily.obs)
            saveRDS(daily.obs, obs.fn) 
        }  # save data in rds 
    }   # end if
  
  })    # end of try()

} # end of subroutine 



# --------------------------------------------------------------------------- #
# inner subroutine to grab observed FLUXNET fluxes, DW, 05/09/2019 
# output fluxes: 30mins, hourly and daily data
# bug fixed when averaging fluxes, remove data whose QC == 3, DW, 05/31/2019
grab.fluxnet.daily.obs <- function(file, site.name, time.zone, min.yr, max.yr, igbp) {

    library(lubridate)

    # ----------------------------------------------------------------------- #
    # read initial data, every 30mins fluxes
    cat('grab.fluxnet.obs(): reading inital 30mins GPP fluxes...it takes time...\n')
    dat <- read.csv(file, sep = ',', header = T)

    # select gap-filled columns
    init.dat <- dat %>% dplyr::select('TIMESTAMP', ends_with('MEAN'), 
                                      starts_with('TA_F'), starts_with('TS_F'),
                                      starts_with('SWC_F'), starts_with('VPD_F'),
                                      starts_with('RH'), ends_with('_MEAN_QC'), 
                                      ends_with('_REF_QC')) 
    
    print(colnames(init.dat))

    cat('REMOVING BAD DATA (QC == 3), QC range before and after: \n')
    print(range(init.dat$NEE_CUT_REF_QC))
    
    # select data within min.yr and max.yr, 
    # remove all poor gapfilled data, QC = 3 and missing data of -9999
    init.dat <- init.dat %>% mutate(timestr = as.character(TIMESTAMP)) %>% 
                             filter(timestr >= paste0(min.yr, '0101'), 
                                    timestr <= paste0(max.yr, '1231')) %>% 
                             filter_at(vars(ends_with('QC')), all_vars(.!= 3)) %>% 
                             filter_at(vars(ends_with('MEAN')), all_vars(.!= -9999)) 

    if (nrow(init.dat) == 0) {
        cat(paste('NO obs data with fine quality (QC != 3) between', min.yr, 
                  'and', max.yr, '\n')); return() }   # end if
    print(range(init.dat$NEE_CUT_REF_QC))
    
    # add time string
    init.dat <- init.dat %>% mutate(site.name = site.name, igbp = igbp, 
                                    LT = as.POSIXct(timestr, format = '%Y%m%d', tz = time.zone), 
                                    DoY = lubridate::yday(LT)) #%>% 

    # ----------------------------------------------------------------------- #
    # now average all columns starting with GPP/NEE/RECO and merge back
    gpp.day.dat  <- init.dat %>% dplyr::select(starts_with('GPP_')) %>% rowMeans()
    reco.day.dat <- init.dat %>% dplyr::select(starts_with('RECO_')) %>% rowMeans()
    nee.day.dat  <- init.dat %>% dplyr::select(starts_with('NEE_')) %>% rowMeans()
    
    avg.dat <- init.dat %>% mutate(GPP = gpp.day.dat, RECO = reco.day.dat, 
                                   NEE = nee.day.dat, site.name = site.name, 
                                   igbp = igbp)

    return(avg.dat) # flux units for daily file is gC/m2/day
}

# end of subroutine
