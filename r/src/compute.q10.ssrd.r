# subroutines to calculate hourly scaling factors, 0-23 UTC
# from hourly ERA5 Tair and SSRD data, DW, 10/24/2019 

#' 12/23/2019, DW, instead of using daily mean ssrd to normalize hourly ssrd, 
#'                 use 4day mean ssrd to match the 4day mean GPP

#' 05/26/2020, DW, need to match radiation data with the exact 4day mean GPP 

compute.q10.ssrd <- function(TA.path, SSRD.path, timestr, site.ext = NULL, 
                             proj.rt = NULL, TA.varname = '2T', 
                             SSRD.varname = c('SSRD', NA)[1], 
                             TA.field = 'ERA5', 
                             SSRD.field = c('ERA5', 'EPIC')[1]) {

    # ------------------------- Prepare incoming SW and TA ------------------- #
    ## Load hourly TA (deg C) and incoming SW from ERA5
    # use nhrs = 23 to grab all hourly variables in a day 
    cat(paste('\ncompute.q10.ssrd(): Loading hourly incoming SW and TA for', timestr, '\n'))
    if (TA.field == 'ERA5') {     # hourly air temp in UTC
        TA.brk <- prep.era5(TA.path, TA.varname, timestr, site.ext, nhrs = 24)

    } else stop(paste('compute.q10.ssrd(): No function available for loading gridded TA from', TA.field))
    # end if TA.field

    # DW, 05/26/2020, before loading SW rad or TA from ERA5
    # create 4day mean interval that matches GPP 
    yr <- substr(timestr, 1, 4)
    date4 <- seq(as.Date(paste0(yr, '-01-01')), as.Date(paste0(yr, '-12-31')), by = '4 days')
    timestr4 <- paste0(format(date4, format = '%Y%m%d'), '00')
    
    # find the start time of the 4-day interval that "timestr" fall into 
    find.timestr <- timestr4[findInterval(timestr, timestr4)]

    # Load hourly surface shortwave radiation downwards (W m-2)
    if (SSRD.field == 'ERA5') {    
        SSRD.brk <- prep.era5(SSRD.path, SSRD.varname, timestr = find.timestr, 
                              site.ext, nhrs = 24 * 4) # now in W/m2

    } else if (SSRD.field == 'EPIC') {
        SSRD.brk <- prep.epic(SSRD.path, timestr = find.timestr, 
                              site.ext, nhrs = 24 * 4)  # initial unit of W m-2

    } else stop(paste('compute.q10.ssrd(): No function available for loading gridded SW downwards from', SSRD.field))
    # end if SW rad.field

    # reproject TA and SWin to 0.05deg to match GPP, RECO
    # default: bilinear interpolation
    TA.pj   <- raster::projectRaster(TA.brk, proj.rt) # in degC
    SSRD.pj <- raster::projectRaster(SSRD.brk, proj.rt)
    SSRD.pj[SSRD.pj < 0] <- 0
    

    # -------------------- Calculate hourly scaling factors ------------------ #
    cat(paste('\ncompute.q10.ssrd(): Calculating Q10, SSRD...\n'))
    # compute Q10 and Tscale from hourly TA (in deg C), from Fisher et al. (2016)
    # in order to downscale daily Reco 
    q10.pj <- 1.5 ^ ((TA.pj - 30) / 10)
    
    # calculate the daily mean of q10 to match daily mean Q10, not daily sum Q10
    # Q10_hr / Q10_daily_avg
    tscale <- q10.pj / mean(q10.pj); names(tscale) <- names(q10.pj)

    # compute Iscale from hourly ERA's SW radiation 
    # DW, 05/26/2020, select the 24 hourly SSRD data from the 4-day data 
    SSRD.dates <- as.POSIXct(names(SSRD.pj), 'UTC', format = 'X%Y.%m.%d.%H.%M.%S')
    SSRD.timestr <- format(SSRD.dates, '%Y%m%d')
    SSRD.sub <- subset(SSRD.pj, which(SSRD.timestr == substr(timestr, 1, 8)))
    
    # calculate the 4-day mean of SSRD, i.e., SSRD_hr / SSRD_daily_avg
    iscale <- SSRD.sub / mean(SSRD.pj); names(iscale) <- names(SSRD.sub)

    return(list(tscale = tscale, iscale = iscale))
}

# end of function



# for sanity check, plot those scaling factors 
if (F) {

    proj.rt[proj.rt == 0] <- NA
    iscale.mk <- raster::mask(iscale, proj.rt); names(iscale.mk) <- seq(0, 23)
    tscale.mk <- raster::mask(tscale, proj.rt); names(tscale.mk) <- seq(0, 23)
    
    i1 <- levelplot(iscale.mk, layout = c(4, 6), at = seq(0, 4.3, 0.1), 
                    xlab = 'LONGITUDE', ylab = 'LATITUDE', max.pixels = 8e6)
    t1 <- levelplot(tscale.mk, layout = c(4, 6), at = seq(0.60, 1.5, 0.05), 
                    xlab = 'LONGITUDE', ylab = 'LATITUDE')
    it <- ggarrange(arrangeGrob(i1), arrangeGrob(t1), ncol = 2, 
                    font.label = list(face = 'plain'),
                    labels = c('a) Iscale for scaling GPP on 1st July 2018', 
                               'b) Tscale for scaling Reco on 1st July 2018'))

    ggsave(it, filename = '../paper3/gmd2020/tscale_iscale_20180701.png',
           width = 12, height = 8)

}









compute.tscale <- function(TA.path, timestr, site.ext = NULL, proj.rt = NULL, 
                           TA.varname = '2T', TA.field = 'ERA5') {

    # ------------------------- Prepare incoming SW and TA ------------------- #
    ## Load hourly TA (deg C) and incoming SW from ERA5
    # use nhrs = 23 to grab all hourly variables in a day 
    cat(paste('\ncompute.q10.ssrd(): Loading hourly incoming SW and TA for', timestr, '\n'))
    if (TA.field == 'ERA5') {     # hourly air temp in UTC
        TA.brk <- prep.era5(TA.path, TA.varname, timestr, site.ext, nhrs = 24)
    } else stop(paste('compute.q10.ssrd(): No function available for loading gridded TA from', TA.field))
    # end if TA.field

    # reproject TA and SWin to 0.05deg to match GPP, RECO
    # default: bilinear interpolation
    TA.pj <- raster::projectRaster(TA.brk, proj.rt)    # in degC

    # -------------------- Calculate hourly scaling factors ------------------ #
    cat(paste('\ncompute.q10.ssrd(): Calculating Q10, SSRD...\n'))
    # compute Q10 and Tscale from hourly TA (in deg C), from Fisher et al. (2016)
    # in order to downscale daily Reco 
    q10.pj <- 1.5 ^ ((TA.pj - 30) / 10)
    
    # calculate the daily mean of q10 to match daily mean Q10, not daily sum Q10
    # Q10_hr / Q10_daily_avg
    tscale <- q10.pj / mean(q10.pj); names(tscale) <- names(q10.pj)
    return(tscale)
}

# end of function


compute.iscale <- function(SSRD.path, timestr, site.ext = NULL, 
                           proj.rt = NULL, SSRD.varname = c('SSRD', NA)[1], 
                           SSRD.field = c('ERA5', 'EPIC')[1]) {

    # ------------------------- Prepare incoming SW ------------------- #
    # Load hourly surface shortwave radiation downwards (W m-2)
    if (SSRD.field == 'ERA5') {    
        SSRD.brk <- prep.era5(SSRD.path, SSRD.varname, timestr, site.ext, nhrs = 24 * 4) # now in W/m2
    } else if (SSRD.field == 'EPIC') {
        SSRD.brk <- prep.epic(SSRD.path, timestr, site.ext, nhrs = 24 * 4)
    } else stop(paste('compute.q10.ssrd(): No function available for loading gridded SW downwards from', SSRD.field))

    # reproject SWin to 0.05deg to match GPP, RECO
    # default: bilinear interpolation
    SSRD.pj <- raster::projectRaster(SSRD.brk, proj.rt); SSRD.pj[SSRD.pj < 0] <- 0

    # -------------------- Calculate hourly scaling factors ------------------ #
    cat(paste('\ncompute.iscale(): Calculating Iscale...\n'))
    # compute Iscale from hourly ERA's SW radiation 
    # calculate the 4-day mean of SSRD, 
    # only return the first 23 layers which correspond to the hourly SSRD for a given day 
    # SSRD_hr / SSRD_daily_avg
    iscale <- SSRD.pj / mean(SSRD.pj); names(iscale) <- names(SSRD.pj)
    iscale.sel <- subset(iscale, 1:24)

    return(iscale.sel)
}

# end of function

