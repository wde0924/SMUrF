#' subroutine to temporally downscale NEE using Tair and SW
#' and spatially downscale using VCF from MOD44B
#' @author Dien Wu, 08/09/2019 

#' 12/23/2019, DW, instead of using daily mean ssrd to normalize hourly ssrd, 
#'                 use 4day mean ssrd to match the 4day mean GPP

downscale.nee.hrly <- function(timestr, gpp.file, reco.file, TA.path, TA.field, 
                               TA.varname, SSRD.path, 
                               SSRD.field = c('ERA5', 'EPIC')[1], 
                               SSRD.varname = c('SSRD', NA)[1]) {

    # ----------------------------- Load GPP and RECO ----------------------- #
    cat(paste('\n\ndownscale.nee.hrly(): Loading daily mean GPP + RECO for', timestr, '\n'))

    # read and clip GPP for every 4 or 8 days
    # select and prepare daily mean GPP in the correct form via crop.smurf.*.r
    gpp.stk <- crop.smurf.gpp(timestr, gpp.file, nhrs = 23, site.ext = NULL)
    reco.stk <- crop.smurf.reco(timestr, reco.file, nhrs = 23, site.ext = NULL)
    mean.gpp.rt <- gpp.stk$GPP_mean; sd.gpp.rt <- gpp.stk$GPP_sd
    mean.reco.rt <- reco.stk$Reco_mean; sd.reco.rt <- reco.stk$Reco_sd

    # site extent will be determined by the overlapped region between GPP and Reco
    # extent(gpp) >= extent(reco.rt), so crop GPP based on RECO, to match RECO
    mean.gpp.int <- raster::intersect(mean.gpp.rt, mean.reco.rt)
    sd.gpp.int <- raster::intersect(sd.gpp.rt, sd.reco.rt)
    site.ext <- extent(mean.gpp.int)
    #rasterOptions(maxmemory = 4E10, memfrac = 0.9)

    # ------------------- Grab hourly Tair and SSRD from ERA5 ---------------- #
    qs.list <- compute.q10.ssrd(TA.path, SSRD.path, timestr, site.ext,
                                proj.rt = mean.reco.rt, TA.varname, 
                                SSRD.varname, TA.field, SSRD.field)
    
    # mean of tscale should be 1
    # mean of 24 hourly iscale != 1, as we used 4-day mean radiation as the normalizer  
    tscale <- qs.list$tscale       
    iscale <- qs.list$iscale      

    # bug fixed, e.g., no ERA5 TA found on 03/11/2018 02UTC, DW, 04/07/2020
    if (nlayers(tscale) != nlayers(iscale)) {
        cat('downscale.nee.hrly(): *** Missing hour in ERA5 file -> dim(t_scale) != dim(i_scale) *** \n')
        
        if (nlayers(tscale) < nlayers(iscale)) 
            iscale <- subset(iscale, which(names(iscale) %in% names(tscale)))
            
        if (nlayers(tscale) > nlayers(iscale)) 
            tscale <- subset(tscale, which(names(tscale) %in% names(iscale)))
    
    }   # end if


    # ------------------ Perform hourly downscaling to fluxes ---------------- #
    cat('downscale.nee.hrly(): Temporally downscaling bio fluxes from daily to hourly...\n')
    mean.reco.stk <- mean.reco.rt * tscale; names(mean.reco.stk) <- names(tscale)
    mean.gpp.stk  <- mean.gpp.int * iscale; names(mean.gpp.stk)  <- names(iscale)
    
    # compute hourly NEE from GPP and Reco 
    mean.nee.stk <- mean.reco.stk - mean.gpp.stk
    names(mean.nee.stk) <- names(mean.reco.stk)

    # return all hourly GPP, Reco and NEE
    hrly.list <- list(mean.gpp.stk, mean.reco.stk, mean.nee.stk)
    names(hrly.list) <- c('hrly_GPP_mean', 'hrly_Reco_mean', 'hrly_NEE_mean')
    
    return(hrly.list)
}   # end of function





# for sanity check, plot those scaling factors 
if (F) {

    selcol <- rev(c('#7F3B08', '#B35806', '#E08214', '#FDB863', '#FEE0B6', '#F7F7F7', 
                    '#D9F0D3', '#A6DBA0', '#5AAE61', '#1B7837', '#00441B'))    
    neeTheme <- rasterTheme(region = selcol)

    plot.nee <- mean.nee.stk; plot.nee[plot.nee == 0] <- NA
    n1 <- levelplot(plot.nee, layout = c(4, 6), at = seq(-55, 55, 5), 
                    par.settings = neeTheme, xlab = 'LONGITUDE', ylab = 'LATITUDE', 
                    max.pixels = 8e6, names.attr = as.character(seq(0, 23)))

    ggsave(arrangeGrob(n1), filename = '../paper3/gmd2020/test_nee_20180701.png',
           width = 6, height = 8)

}


