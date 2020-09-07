#' a funcion inside predReco() to predict Reco based on pre-trained NN
#' @author DW, 08/01/2019

#' call several inner functions (nerualnet.reco.urban and nerualnet.reco)
#' see those functions in Rscript neuralnet_functions.R

#' update by DW:
#' 08/02/2019, add biomes-specific Reco uncertainty
#' 09/10/2019, use updated Land Cover tif generated from main_script_GPP.r
#' 12/02/2019, for a few grid cells near water bodies whose temps are NA, 
#'             use a simple linear relation between GPP and Reco, ~ about half
#' 03/02/2020, instead of gap-filled land cover, reuse vegetated fractions. 
#' 04/20/2020, add C3C4 partition and combine OSHR and CSHR as SHR
#' 04/30/2020, add Daymet + NLDAS for predicting Reco for CONUS

# ---------------------------------------------------------------------------- #
predReco_loop_neuralnet <- function(prep.list, nn.dir, 
                                    nn.pattern = c('daymet_nldas', 'era5')[1]) {
    
    input.brk <- prep.list$input.brk 
    lc.rt     <- prep.list$lc.rt
    fveg.stk  <- prep.list$fveg.stk
    c3c4.stk  <- prep.list$c3c4.stk 

    # search for all available NN models and get CV for Reco error
    nn.fns <- list.files(nn.dir, nn.pattern, full.names = T)
    if (nn.pattern == 'daymet_nldas') 
        nn.fns <- list.files(nn.dir, pattern = 'daymet\\+nldas', full.names = T)


    ## get land cover # with their abbreviation ------------------------------- #
    # in the order of ENF, EBF, DNF, DBF, MF, CSHR, OSHR, WSAV, SAV, GRA, 
    #                 WET, CRO, URB, CRO/NVM, SI, BAR, WAT, UNC, WGT, CRO_C3, 
    #                 CRO_C4, SHR
    cv.era5 = c(0.419, 0.381, 0.726, 0.421, 0.462, NA, NA, 0.34, 0.352, 0.445, 
                0.412, 0.458,    NA,    NA,    NA, NA, NA,   NA,    NA, 0.472, 
                0.390, 0.559)
    
    # No EBF, DNF sites for US when using daymet+NLDAS, so CV is NA
    cv.daymet.nldas <- c(0.384,    NA, NA, 0.352, 0.65, NA, 0.613, 0.531, NA, 0.533, 
                          0.34, 0.428, NA, 0.428,   NA, NA,    NA,    NA, NA, 0.497, 
                          0.37,    NA)    
    if (nn.pattern == 'era5') cv.vec = cv.era5 
    if (nn.pattern == 'daymet_nldas') cv.vec = cv.daymet.nldas 
    igbp.info <- get.igbp.col() %>% dplyr::select(c('val', 'name')) %>%
                 tibble::add_row(val = 23, name = 'SHR') %>% mutate(cv = cv.vec)


    ## modify lc.rt accordingly ----------------------------------------------
    lc.new <- lc.rt 

    # FOR the GLOBE: remove CRO or CRO/NVM, use CRO_C3 and CRO_C4 instead
    # combine open and closed SHR as SHR due to low # of flux sites available
    lc.new[lc.new == 14] <- 12                # let CRO/NVM be CRO (# = 12)
    if (nn.pattern == 'era5') lc.new[lc.new == 6 | lc.new == 7] <- 23  
    
    # FOR CONUS: let CSHR be OSHR and SAV be WSAV when using daymet+nldas
    if (nn.pattern == 'daymet_nldas') {
        lc.new[lc.new == 6] <- 7    # 7 for OSHR 
        lc.new[lc.new == 9] <- 8    # 8 for WSAV
    }   # end if

    # ------------------------------------------------------------------------ #
    # normalize TA, TS, and GPP, 
    # please do not change these if using default pretrained NN models
    min.df <- data.frame('Reco' = -1, 'TA' = -35, 'TS' = -30, 'GPP' = -1)
    max.df <- data.frame('Reco' = 20, 'TA' =  40, 'TS' =  40, 'GPP' = 30)

    # calculate the coeffient to convert normalized scales back to regular scales
    # e.g., x = x_norm * a + b where x_norm fall within [0, 1]
    a <- max.df - min.df; b <- min.df 

    # normalize based on prescribed min-max boundary --------------------------
    cat('predReco_loop_neuralnet(): normalizing input variables...\n')
    norm.brk <- input.brk
    norm.brk$TA  <- (norm.brk$TA - b$TA) / a$TA
    norm.brk$TS  <- (norm.brk$TS - b$TS) / a$TS
    norm.brk$GPP <- (norm.brk$GPP - b$GPP) / a$GPP

    # if use NN model that trained based on all biomes
    #nn.model <- readRDS(nn.fn[grepl('ALL', nn.fn)])
    #Reco <- raster::predict(norm.brk, nn.model$nn)
    #Reco <- Reco * a$Reco + b$Reco
    
    # ------------------------------------------------------------------------ #
    # biomes that will get a Reco value
    igbp.abbr <- (igbp.info %>% filter(val %in% sort(unique(lc.new))))$name
    
    # skip biomes of EBF and DNF IF using daymet and NLDAS
    # Reco for grid of EBF and DNF will be NA
    if (nn.pattern == 'daymet_nldas') 
        igbp.abbr <- igbp.abbr[!igbp.abbr %in% c('EBF', 'DNF')]


    # ------------------------------------------------------------------------ #
    for (i in 1 : length(igbp.abbr)) {
        
        tmp.name <- igbp.abbr[i]
        tmp.info <- igbp.info[igbp.info$name == tmp.name, ]
        cat(paste0('predReco_loop_neuralnet(): working on ', tmp.name, ';'))
        
        # grab IGBP info and select land cover map for the desired land types
        lc.tmp <- lc.new; lc.tmp[!lc.tmp %in% tmp.info$val] <- NA

        # ------------------------------------------------------------------- #
        # force Reco and uncertainties for 500 m grids of SI, BAR, WAT, UNC as zero
        if (tmp.name %in% c('WAT', 'SI', 'BAR', 'UNC')) {
            tmp.reco.mean <- lc.tmp 
            tmp.reco.mean[!is.na(tmp.reco.mean)] <- 0
            tmp.reco.sd <- tmp.reco.mean

        } else if (tmp.name == 'URB') {
            
            # select input variables for a certain land cover type
            x <- raster::mask(norm.brk, lc.tmp)
            y.urb <- neuralnet.reco.urban(x, igbp.info, fveg.stk, nn.fns, nn.pattern)

            # convert Reco back to regular scale
            tmp.reco.mean <- y.urb$y * a$Reco + b$Reco
            tmp.reco.sd   <- tmp.reco.mean * y.urb$y.cv 

        } else if (tmp.name == 'CRO') {
            
            # select input variables for a certain land cover type
            x <- raster::mask(norm.brk, lc.tmp)

            cat('Predicting Reco for CROPs assuming C3 or C4 crops...\n')
            y.cro3 <- neuralnet.reco(x, nn.fns[grepl('_CRO_C3', nn.fns)])
            y.cro4 <- neuralnet.reco(x, nn.fns[grepl('_CRO_C4', nn.fns)])
            y <- y.cro3 * c3c4.stk$C3_ratio + y.cro4 * c3c4.stk$C4_ratio
            
            cv.cro3 <- igbp.info[igbp.info$name == 'CRO_C3', 'cv']
            cv.cro4 <- igbp.info[igbp.info$name == 'CRO_C4', 'cv']
            y.cv <- sqrt((cv.cro3 * c3c4.stk$C3_ratio) ^ 2 + (cv.cro4 * c3c4.stk$C4_ratio) ^ 2)
            
            # convert Reco back to regular scale
            tmp.reco.mean <- y * a$Reco + b$Reco
            tmp.reco.sd   <- tmp.reco.mean * y.cv

        } else {

            # select input variables for a certain land cover type
            x <- raster::mask(norm.brk, lc.tmp)

            # make prediction and convert Reco back to regular scale
            y <- neuralnet.reco(x, nn.fns[grepl(paste0('_', tmp.name), nn.fns)])
            y.cv <- tmp.info$cv 
            tmp.reco.mean <- y * a$Reco + b$Reco
            tmp.reco.sd   <- tmp.reco.mean * y.cv

        }   # end if

        cat(paste0(' min/max Reco at 500 m: ', signif(minValue(tmp.reco.mean), 3), 
                   '/', signif(maxValue(tmp.reco.mean), 3), '\n'))

        # merge Reco from different biomes
        if (i == 1) { 
            reco.mean <- tmp.reco.mean; reco.sd <- tmp.reco.sd
        } else {
            reco.mean <- raster::merge(reco.mean, tmp.reco.mean)
            reco.sd   <- raster::merge(reco.sd,   tmp.reco.sd)
        }   # end if

        gc() 
    }   # end for i 


    # for water, ice, and barren land, force Reco as Zero instead of NA 
    # this will avoid huge Reco number over grid cells next to water
    pred.stk <- stack(reco.mean, reco.sd)
    #pred.stk[is.na(pred.stk)] <- 0
    names(pred.stk) <- c('reco.mean', 'reco.sd')

    cat('predReco_loop_neuralnet(): END of this function of predicting Reco at 500m')
    return(pred.stk)
}   
# end of predReco.loop()

