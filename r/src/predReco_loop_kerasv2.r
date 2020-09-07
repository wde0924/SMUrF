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
predReco_loop_kerasv2 <- function(timestr, prep.list, nn.file) {
    
    input.brk <- prep.list$input.brk 
    lc.rt     <- prep.list$lc.rt
    fveg.stk  <- prep.list$fveg.stk
    c3c4.stk  <- prep.list$c3c4.stk 

    # get land cover # with their abbreviation ------------------------------- #
    # remove CRO or CRO/NVM, use CRO_C3 and CRO_C4 instead
    # combine open and closed SHR as SHR due to low # of  flux sites available
    # i.e., 1 ENF; 2 EBF; 3 DNF; 4 DBF; 5 MF; 6 CSHR; 7 OSHR; 8 WSAV; 9 SAV; 
    #       10 GRA; 11 WET; 12 CRO; 13 URB; 21 CRO_C3; 22 CRO_C4
    igbp.info <- get.igbp.col() %>% dplyr::select(c('val', 'name')) 

    # modify lc.rt, let CRO/NVM be CRO (# = 12)
    lc.new <- lc.rt; lc.new[lc.new == 14] <- 12
    #lc.new[lc.new %in% c(15, 16, 17, 254, 20)] <- NA
    names(lc.new) <- 'LC'

    # ------------------------------------------------------------------------ #
    # normalize TA, TS, and GPP, 
    # please do not change these values if using default pretrained NN models
    # ------------------------------------------------------------------------ #   
    min.df <- data.frame('Reco' =  0, 'TA' = -52, 'TS' = -47, 'GPP' = -2)
    max.df <- data.frame('Reco' = 20, 'TA' =  38, 'TS' =  41, 'GPP' = 19)

    # calculate the coeffient to convert normalized scales back to regular scales
    # e.g., x = x_norm * a + b where x_norm fall within [0, 1]
    a <- max.df - min.df; b <- min.df 

    # normalize based on prescribed min-max boundary --------------------------
    cat('predReco_loop_kerasv2(): normalizing input variables...\n')
    norm.stk <- stack(input.brk, lc.new)
    norm.stk$TA  <- (norm.stk$TA - b$TA) / a$TA
    norm.stk$TS  <- (norm.stk$TS - b$TS) / a$TS
    norm.stk$GPP <- (norm.stk$GPP - b$GPP) / a$GPP
    
    # convert to data frame and matrix 
    norm.df <- as.data.frame(norm.stk, xy = T) %>% 
               left_join(igbp.info, by = c('LC' = 'val')) %>% mutate(Reco = 0)
    norm.df$row.indx <- 1 : nrow(norm.df)


    # ------------------------------------------------------------------------ #     
    # add month and season indics, if southern hemisphere, flip month
    # ------------------------------------------------------------------------ #   
    cat('predReco_loop_kerasv2(): converting month/season to one-hot encoding...\n')
    mon_init = as.numeric(substr(timestr, 5, 6))
    mon_flip = ifelse(min(norm.df$y) < 0, mon_init - 6, mon_init)
    mon_flip = ifelse(min(norm.df$y) < 0 & mon_flip <= 0, mon_flip + 12, mon_flip)
    
    # 1 spring, 2 summer, 3 fall, 4 winter
    season_indx = ifelse(mon_flip %in% c(3, 4, 5), 1,  
                  ifelse(mon_flip %in% c(6, 7, 8), 2,  
                  ifelse(mon_flip %in% c(9, 10, 11), 3, 4)))
 
    # convert month and season indices to one-hot encoding
    pred.df <- norm.df %>% dplyr::select('x', 'y', 'GPP', 'TA', 'TS', 'LC', 'row.indx') %>% 
               mutate(mon_1 = 0, mon_2 = 0, mon_3 = 0, mon_4 = 0, mon_5 = 0, 
                      mon_6 = 0, mon_7 = 0, mon_8 = 0, mon_9 = 0, mon_10 = 0, 
                      mon_11 = 0, mon_12 = 0, season_1 = 0, season_2 = 0, 
                      season_3 = 0, season_4 = 0)
    
    pred.df[, colnames(pred.df) == paste0('mon_', mon_flip)] <- 1
    pred.df[, colnames(pred.df) == paste0('season_', season_indx)] <- 1
    pred.df$Reco <- 0

    # load model 
    cat('predReco_loop_kerasv2(): loading NN model...\n')
    model <- load_model_hdf5(nn.file)


    # ------------------------------------------------------------------------ #     
    # STEP 1) select and work on natural biomes except for CROP and URBAN
    # ------------------------------------------------------------------------ # 
    cat('\n\npredReco_loop_kerasv2(): # ------- STEP 2.1: working on 11 natural biomes ------- #\n')
    natural.df <- pred.df %>% filter(LC <= 11, !is.na(GPP))
    natural.mtrx <- as.matrix(natural.df)

    # convert Land cover type to one-hot encoding, 0 to num.classes
    cat('Land use type one-hot encoding...\n')
    lc.mtrx <- keras::to_categorical(y = natural.df$LC - 1, num_classes = 11)
    colnames(lc.mtrx) <- paste0('igbp_', seq(1, 11))
    lc.mtrx <- cbind(lc.mtrx, igbp_21 = 0, igbp_22 = 0)
    natural.mtrx <- cbind(natural.mtrx, lc.mtrx)
    natural.mtrx <- natural.mtrx[, which(!colnames(natural.mtrx) %in% c('LC', 'row.indx', 'x', 'y', 'Reco'))]

    # make prediction and put predicted Reco back to pred.df 
    if (nrow(natural.mtrx) == 0) stop('predReco_loop_kerasv2(): empty matrix for Reco prediction, please check...\n')
    cat('Predicting 11 natural biomes, it usually takes ~15 mins given large amount of data...\n')
    natural.reco <- as.numeric(model %>% predict(natural.mtrx))
    pred.df[which(pred.df$row.indx %in% natural.df$row.indx), 'Reco'] <- natural.reco


    # ------------------------------------------------------------------------ # 
    # STEP 2) select and work on CROP
    # ------------------------------------------------------------------------ # 
    cat('\n\npredReco_loop_kerasv2(): # ------- STEP 2.2: working on croplands ------- #\n')
    crop.df <- pred.df %>% filter(LC == 12, !is.na(GPP)) %>% 
               mutate(igbp_1 = 0, igbp_2 = 0, igbp_3 = 0, igbp_4 = 0, igbp_5 = 0,
                      igbp_6 = 0, igbp_7 = 0, igbp_8 = 0, igbp_9 = 0, igbp_10 = 0, 
                      igbp_11 = 0, igbp_21 = 0, igbp_22 = 0)
    crop.mtrx <- as.matrix(crop.df %>% dplyr::select(-c('LC', 'row.indx', 'Reco', 'x', 'y')))
    cro3.mtrx = crop.mtrx; cro3.mtrx[, 'igbp_21'] <- 1
    cro4.mtrx = crop.mtrx; cro4.mtrx[, 'igbp_22'] <- 1

    # calculate the weighted mean crop reco based on C3 C4 ratios
    cat('Predicting Reco for CROPs assuming C3 crops...\n')
    cro3.reco <- as.numeric(model %>% predict(cro3.mtrx))

    cat('Predicting Reco for CROPs assuming C4 crops...\n')
    cro4.reco <- as.numeric(model %>% predict(cro4.mtrx))

    # interpolate C4 ratio onto 500m CRO
    cat('Calculating weighted mean Reco for CROPs...\n')
    c4.frac <- raster::extract(c3c4.stk$C4_ratio, crop.df[, c('x', 'y')])
    c4.frac[is.na(c4.frac)] <- 0
    
    # calculate weighted mean Reco for crops
    crop.reco <- cro3.reco * (1 - c4.frac) + cro4.reco * c4.frac
    pred.df[which(pred.df$row.indx %in% crop.df$row.indx), 'Reco'] <- crop.reco


    # ------------------------------------------------------------------------ # 
    # STEP 3) select and work on URBAN grids
    # ------------------------------------------------------------------------ # 
    cat('\n\npredReco_loop_kerasv2(): # ------- STEP 2.3: working on URBAN ------- #\n')
    urban.df <- pred.df %>% filter(LC == 13, !is.na(GPP)) %>% 
                mutate(igbp_1 = 0, igbp_2 = 0, igbp_3 = 0, igbp_4 = 0, igbp_5 = 0,
                       igbp_6 = 0, igbp_7 = 0, igbp_8 = 0, igbp_9 = 0, igbp_10 = 0, 
                       igbp_11 = 0, igbp_21 = 0, igbp_22 = 0)
                
    urban.mtrx <- as.matrix(urban.df %>% dplyr::select(-c('LC', 'row.indx', 'Reco', 'x', 'y')))

    # interpolate veg fractions at urban grid 
    cat('Extracting vegetation fractions for URBAN areas...\n')
    fveg.mtrx <- raster::extract(fveg.stk, urban.df[, c('x', 'y')])

    cat('Predicting Reco for URBAN as various veg types...\n')
    enf.mtrx = urban.mtrx; enf.mtrx[, 'igbp_1'] <- 1
    ebf.mtrx = urban.mtrx; ebf.mtrx[, 'igbp_2'] <- 1
    dnf.mtrx = urban.mtrx; dnf.mtrx[, 'igbp_3'] <- 1
    dbf.mtrx = urban.mtrx; dbf.mtrx[, 'igbp_4'] <- 1
    mf.mtrx  = urban.mtrx; mf.mtrx[, 'igbp_5'] <- 1
    shr.mtrx = urban.mtrx; shr.mtrx[, 'igbp_7'] <- 1
    gra.mtrx = urban.mtrx; gra.mtrx[, 'igbp_10'] <- 1

    # ------------------------------------------------------------------------ #  
    enf.reco <- as.numeric(model %>% predict(enf.mtrx))
    ebf.reco <- as.numeric(model %>% predict(ebf.mtrx))
    dnf.reco <- as.numeric(model %>% predict(dnf.mtrx))
    dbf.reco <- as.numeric(model %>% predict(dbf.mtrx))
    mf.reco  <- as.numeric(model %>% predict(mf.mtrx))
    shr.reco <- as.numeric(model %>% predict(shr.mtrx))
    gra.reco <- as.numeric(model %>% predict(gra.mtrx))

    # weight Reco based on fveg.stk being proccessed during GPP estimates
    cat('Calculating weighted mean Reco for URBAN areas...\n')
    urban.reco <- gra.reco * fveg.mtrx[, 'GRA'] + shr.reco * fveg.mtrx[, 'OSHR'] + 
                  dbf.reco * fveg.mtrx[, 'DBF'] + dnf.reco * fveg.mtrx[, 'DNF'] + 
                  ebf.reco * fveg.mtrx[, 'EBF'] + enf.reco * fveg.mtrx[, 'ENF'] + 
                  mf.reco * fveg.mtrx[, 'MF']
    pred.df[which(pred.df$row.indx %in% urban.df$row.indx), 'Reco'] <- urban.reco


    # ------------------------------------------------------------------------ #
    # STEP 4) transform back to regular scales
    # ------------------------------------------------------------------------ #
    #norm.df[which(norm.df$row.indx %in% pred.df$row.indx), 'Reco'] <- pred.df$Reco 
    cat('\n\npredReco_loop_kerasv2(): # ------- STEP 2.4: Converting Reco back to regular scales -------- #\n')
    output.df <- pred.df %>% 
                 # if TA or TS is NA, use half of GPP as Reco 
                 mutate(#Reco = ifelse(is.na(TA) | is.na(TS), 0.5 * GPP, Reco), 

                 # convert GPP and Reco back to regular scale
                        GPP = GPP * a$GPP + b$GPP, Reco = Reco * a$Reco + b$Reco, 
                        Reco_sd = Reco * 0.455) 

    # reform dataframe back to rasterLayer
    output.brk <- rasterFromXYZ(output.df[, c('x', 'y', 'Reco', 'Reco_sd')])
    crs(output.brk) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    names(output.brk) <- list('reco.mean', 'reco.sd')
    
    gc() 
    cat('predReco_loop_kerasv2(): END of this function of predicting Reco at 500m')
    return(output.brk)
}   

# end of function

