# function to prepare input variables for Reco prediction, DW, 03/04/2020
# no need to go down to 500 m if we did not differentiate biomes 
# no need to grab vegetated fractions
# only downscale to 0.05 degree is enough

prep.input4reco.keras <- function(reg.ext, reg.path, reg.name, timestr, 
                                  TA.field, TA.path, TA.varname,
                                  TS.field, TS.path, TS.varname, gpp.file) {

    # ---------------------- Step 1.1 PREPARE Tair --------------------------- #
    ## numbers of hrs (forward or backward) for aggregating and averaging from 00UTC
    # air and soil temperatures (to estimate Reco fluxes)
    # if nhrs is NULL, calculate the Reco right at 00UTC 
    nhrs <- 23           # simply use daily mean Tair, + 23 hrs
    yr <- substr(timestr, 1, 4)

    ## calculate n-day daily mean Tair (deg C) from daymet via prep.tair.daymet()
    # output is a rasterlayer with approximate grid spacing
    if (TA.field == 'daymet') {
        TA.brk <- prep.tair.daymet(TA.path, TA.varname, timestr, reg.ext, nhrs)
    } else if (TA.field == 'ERA5') {
        TA.brk <- prep.era5(TA.path, TA.varname, timestr, reg.ext, nhrs)
    } else {
        stop(paste('prep.input4reco(): No function available for loading gridded Tair for', 
                    TA.field))
    }   # end if TA.field
    
    # take the average if more than 1 layer exists
    if (raster::nlayers(TA.brk) > 1) { mean.TA.rt <- mean(TA.brk) 
    } else mean.TA.rt <- TA.brk     # already in degC for Tair


    # ---------------------- Step 1.2 PREPARE Tsoil -------------------------- #
    ## get hourly Tsoil data by calling prep.tsoil.*(), add ERA5 Tsoil, DW, 07/01/2019
    # then calculate mean temperature
    if (TS.field == 'NLDAS') {
        TS.brk <- prep.soil.nldas(TS.path, TS.varname, timestr, reg.ext, nhrs)
    } else if (TS.field == 'ERA5') {
        TS.brk <- prep.era5(TS.path, TS.varname, timestr, reg.ext, nhrs)
    } else {
        stop(paste('prep.input4reco(): No function available for loading gridded Tsoil for', 
                    TS.field))
    }   # end if TS.field
       
    if (raster::nlayers(TS.brk) > 1) { mean.TS.rt  <- mean(TS.brk)
    } else mean.TS.rt <- TS.brk


    # ---------------------- Step 1.3 PREPARE GPP -------------------------- #
    # read and clip GPP for every 4 or 8 days
    gpp.file <- gpp.file[grepl(yr, gpp.file)]
    if (length(gpp.file) == 0) 
        stop('prep.input4reco(): NO SIF-based GPP found, see `main_script_GPP.r`\n')
  
    # select and prepare daily mean GPP in the correct form via prep.gpp()
    gpp.stk <- crop.smurf.gpp(timestr, gpp.file, nhrs, reg.ext)
    #rasterOptions(maxmemory = 4E10, memfrac = 0.9)

    # ---------------------- Step 1.5 REPROJECTION -------------------------- #
    ## re-project mean Tsoil and Tair to 500m land cover data with bilinear interpolation 
    cat('prep.input4reco(): Reprojecting Tair, Tsoil, GPP to 500m to match MCD IGBP...\n')
    TA.pj <- raster::projectRaster(mean.TA.rt, gpp.stk); gc()
    TS.pj <- raster::projectRaster(mean.TS.rt, gpp.stk); gc()

    # merge all input as a rasterBrick interpolated at 500m 
    input.brk <- brick(TA.pj, TS.pj, gpp.stk$GPP_mean)
    names(input.brk) <- c('TA', 'TS', 'GPP')

    return(input.brk)

}
