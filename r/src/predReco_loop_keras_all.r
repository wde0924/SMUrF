#' a funcion inside predReco() to predict Reco based on pre-trained NN
#' @author DW, 08/01/2019

#' add biomes-specific Reco uncertainty, DW, 08/02/2019

#' update by DW:
#' 09/10/2019, use updated Land Cover tif generated from main_script_GPP.r
#' 09/19/2019, increase the Tair for urban by 2 deg C. 
#' 12/02/2019, for a few grid cells near water bodies whose temps are NA, 
#'             use a simple linear relation between GPP and Reco, ~ about half
#' 03/02/2020, instead of gap-filled land cover, reuse vegetated fractions. 
#' 04/02/2020, treat all biomes together, no need for MODIS IGBP
#' 04/06/2020 need to force Reco as NA for water, ice, and barren areas


# ---------------------------------------------------------------------------- #
predReco_loop_keras_all <- function(input.brk, TA.field, TS.field, lc.rt, smurf_wd) {
    
    # search for all available NN models and get CV for Reco error
    nn.files <- list.files(file.path(smurf_wd, 'data'), '.h5', full.names = T)

    # normalize TA, TS, and GPP, 
    # please do not change these if using default pretrained NN models
    min.df <- data.frame('Reco' = 0, 'TA' = -35, 'TS' = -30, 'GPP' = -1)
    max.df <- data.frame('Reco' = 20, 'TA' = 40, 'TS' = 40, 'GPP' = 30)

    # calculate the coeffient to convert normalized scales back to regular scales
    # e.g., x = x_norm * a + b where x_norm fall within [0, 1]
    a <- max.df - min.df; b <- min.df 

    # normalize based on prescribed min-max boundary --------------------------
    input.brk$TA  <- (input.brk$TA - b$TA) / a$TA
    input.brk$TS  <- (input.brk$TS - b$TS) / a$TS
    input.brk$GPP <- (input.brk$GPP - b$GPP) / a$GPP
    input.df   <- as.data.frame(input.brk, xy = T)
    input.mtrx <- as.matrix(input.df %>% dplyr::select(-c('x', 'y')))
    
    # read h5 file for NN models and predict Reco for natural biomes
    if (TA.field == 'ERA5' & TS.field == 'ERA5') {
        cat('predReco_loop_keras(): use the NN model trained based on ERA5...\n')
        model <- keras::load_model_hdf5(nn.files[grepl('era5', nn.files)])
        cv.reco <- 0.452

    } else if (TA.field == 'daymet' & TS.field == 'NLDAS') {
        cat('predReco_loop_keras(): use the NN model trained based on daymet + NLDAS...\n')
        model <- keras::load_model_hdf5(nn.files[grepl('NLDAS', nn.files)])
        cv.reco <- 0.501

    } else {
        cat('predReco_loop_keras(): use the NN model trained based on FLUXNET...\n')
        model <- keras::load_model_hdf5(nn.files[grepl('FLUXNET', nn.files)])
        cv.reco <- 0.426
    }   # end if load model 

    # predict Reco 
    input.df$Reco <- as.numeric(model %>% predict(input.mtrx))
    
    # transform back to regular scales and reform dataframe back to rasterLayer
    output.df <- input.df %>% mutate(TA = TA * a$TA + b$TA, 
                                     TS = TS * a$TS + b$TS, 
                                     GPP = GPP * a$GPP + b$GPP, 
                                     Reco = Reco * a$Reco + b$Reco, 
                                     Reco_sd = Reco * cv.reco)    
                              # fractional uncertainty of 49.1%

    output.brk <- rasterFromXYZ(output.df)
    crs(output.brk) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    
    # interpolate Reco at 500m ------------------------------------------------
    # force Reco over grid cells with water, barren, and ice as zero 
    cat('predReco_loop_keras(): forcing Reco for water, barren, and ice as NA...\n')

    # remove biomes e.g., 15, 16, 17 for SI, BAR, WAT, 254 for UNC
    rm.rt <- lc.rt; rm.rt[rm.rt %in% c(15, 16, 17, 254)] <- NA

    reco.stk <- output.brk[[which(grepl('Reco', names(output.brk)))]]
    reco.ds  <- projectRaster(reco.stk, lc.rt)  # downscale to 500m first
    reco.rm  <- raster::mask(reco.ds, rm.rt)    # let Reco be NA for water, ice...
    reco.agg <- aggregate(reco.rm, fact = 12, na.rm = T)    # aggregate to 0.05deg

    return(reco.agg)
}   

# end of function

