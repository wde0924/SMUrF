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
predReco_loop_keras <- function(prep.list, nn.dir, nn.pattern) {
    
    input.brk <- prep.list$input.brk 
    lc.rt     <- prep.list$lc.rt
    fveg.stk  <- prep.list$fveg.stk
    c3c4.stk  <- prep.list$c3c4.stk 

    # search for all available NN models and get CV for Reco error
    nn.file <- list.files(nn.dir, nn.pattern, full.names = T)

    # get land cover # with their abbreviation ------------------------------- #
    # remove CRO or CRO/NVM, use CRO_C3 and CRO_C4 instead
    # combine open and closed SHR as SHR due to low # of  flux sites available
    igbp.info <- get.igbp.col() %>% dplyr::select(c('val', 'name')) %>%
                 mutate(cv = c(0.419, 0.381, 0.726,  0.421,  0.462, 
                               0.559, 0.559,  0.34,  0.352,  0.445, 
                               0.412, 0.458,    NA,     NA,     NA, 
                                  NA,    NA,    NA,     NA,     NA, 
                                  NA))        # CV using ERA5

    # modify lc.rt, let CRO/NVM be CRO (# = 12)
    lc.new <- lc.rt; lc.new[lc.new == 14] <- 12; names(lc.new) <- 'LC'
    
    # ------------------------------------------------------------------------ #
    # normalize TA, TS, and GPP, 
    # please do not change these if using default pretrained NN models
    min.df <- data.frame('Reco' = -1, 'TA' = -20, 'TS' = -20, 'GPP' = -1)
    max.df <- data.frame('Reco' = 20, 'TA' = 40, 'TS' = 40, 'GPP' = 30)

    # calculate the coeffient to convert normalized scales back to regular scales
    # e.g., x = x_norm * a + b where x_norm fall within [0, 1]
    a <- max.df - min.df; b <- min.df 

    # normalize based on prescribed min-max boundary --------------------------
    cat('predReco_loop_keras(): normalizing input variables...\n')
    norm.stk <- stack(input.brk, lc.new)
    norm.stk$TA  <- (norm.stk$TA - b$TA) / a$TA
    norm.stk$TS  <- (norm.stk$TS - b$TS) / a$TS
    norm.stk$GPP <- (norm.stk$GPP - b$GPP) / a$GPP
    
    # convert to data frame and matrix 
    norm.df <- as.data.frame(norm.stk, xy = T) %>% 
               left_join(igbp.info, by = c('LC' = 'val')) %>% mutate(Reco = 0)
    norm.df$row.indx <- 1 : nrow(norm.df)

    # biomes that will get a Reco value, SI, BAR, WAT, UNC will be neglected ---
    igbp.abbr <- c('ENF', 'EBF', 'DNF', 'DBF', 'MF', 'CSHR', 'OSHR', 'WSAV', 
                   'SAV', 'GRA', 'WET', 'CRO', 'URB')

    # ------------------------------------------------------------------------ #     
    for (i in 1 : length(igbp.abbr)) {

        tmp.name <- igbp.abbr[i]
        tmp.info <- igbp.info[igbp.info$name == tmp.name, ]

        # select input variables for a certain land cover type
        x.df <- norm.df %>% filter(name == tmp.name) 
        if (nrow(x.df) == 0) next 
        cat(paste('predReco_loop_keras(): working on', tmp.name,
                   signif(nrow(x.df) / nrow(norm.df), 3) * 100, '% of the data\n'))
        
        x.mtrx <- as.matrix(x.df %>% dplyr::select(c('TA', 'TS', 'GPP')))

        # load model and make prediction
        pred.inner <- function(veg, x.mtrx, nn.file) {
            model <- keras::load_model_hdf5(nn.file[grepl(paste0('_', veg), nn.file)])
            as.numeric(model %>% predict(x.mtrx))
        }

        if (tmp.name == 'URB') {
            
            cat('Urban Reco estimates take a while...\n')
            y.gra <- pred.inner('GRA', x.mtrx, nn.file)
            y.shr <- pred.inner('OSHR', x.mtrx, nn.file)
            y.dbf <- pred.inner('DBF', x.mtrx, nn.file)
            y.dnf <- pred.inner('DNF', x.mtrx, nn.file)
            y.ebf <- pred.inner('EBF', x.mtrx, nn.file)
            y.enf <- pred.inner('ENF', x.mtrx, nn.file)
            y.mf  <- pred.inner('MF', x.mtrx, nn.file); gc() 

            # get fractional uncertainty for GRA, DBF, EBF, and ENF 
            gra.cv <- igbp.info[igbp.info$name == 'GRA', 'cv']
            shr.cv <- igbp.info[igbp.info$name == 'OSHR', 'cv']
            dbf.cv <- igbp.info[igbp.info$name == 'DBF', 'cv']
            dnf.cv <- igbp.info[igbp.info$name == 'DNF', 'cv']
            ebf.cv <- igbp.info[igbp.info$name == 'EBF', 'cv']
            enf.cv <- igbp.info[igbp.info$name == 'ENF', 'cv']
            mf.cv  <- igbp.info[igbp.info$name == 'MF', 'cv']

            # weight Reco based on fveg.stk being proccessed during GPP estimates
            # interpolate veg fractions at urban grid 
            fveg.mtrx <- raster::extract(fveg.stk, x.df[, c('x', 'y')])
            y.reco <- y.gra * fveg.mtrx[, 'GRA'] + y.shr * fveg.mtrx[, 'OSHR'] + 
                      y.dbf * fveg.mtrx[, 'DBF'] + y.dnf * fveg.mtrx[, 'DNF'] + 
                      y.ebf * fveg.mtrx[, 'EBF'] + y.enf * fveg.mtrx[, 'ENF'] + 
                      y.mf * fveg.mtrx[, 'MF']

            y.cv <- sqrt((gra.cv * fveg.mtrx[, 'GRA']) ^ 2 + 
                         (shr.cv * fveg.mtrx[, 'OSHR']) ^ 2 + 
                         (dbf.cv * fveg.mtrx[, 'DBF']) ^ 2 + 
                         (dnf.cv * fveg.mtrx[, 'DNF']) ^ 2 + 
                         (ebf.cv * fveg.mtrx[, 'EBF']) ^ 2 + 
                         (enf.cv * fveg.mtrx[, 'ENF']) ^ 2 + 
                         (mf.cv * fveg.mtrx[, 'MF']) ^ 2) 

        } else if (tmp.name == 'CRO') {
            
            # calculate the weighted mean crop reco based on C3 C4 ratios
            cat('performing C3-C4 partitioning which takes a while...\n')
            y.cro3 <- pred.inner('CRO_C3', x.mtrx, nn.file)
            y.cro4 <- pred.inner('CRO_C4', x.mtrx, nn.file)

            # interpolate C4 ratio onto 500m CRO
            c4.frac <- raster::extract(c3c4.stk$C4_ratio, x.df[, c('x', 'y')])
            c4.frac[is.na(c4.frac)] <- 0
            y.reco <- y.cro3 * (1 - c4.frac) + y.cro4 * c4.frac
            y.cv <- tmp.info$cv 

        } else {

            # load model and make prediction
            y.reco <- pred.inner(tmp.name, x.mtrx, nn.file)
            y.cv <- tmp.info$cv 

        }   # end if tmp.name

        # put Reco and fractional uncertainties back to norm.df
        norm.df[which(norm.df$row.indx %in% x.df$row.indx), 'Reco'] <- y.reco
        norm.df[which(norm.df$row.indx %in% x.df$row.indx), 'cv'] <- y.cv

        print(range(x.df$GPP * 31 - 1, na.rm = T))
        print(range(y.reco * 21 - 1, na.rm = T))

    }   # end for i

    
    # ------------------------------------------------------------------------ #
    # transform back to regular scales
    output.df <- norm.df %>% mutate(TA = TA * a$TA + b$TA, 
                                    TS = TS * a$TS + b$TS, 
                                    GPP = GPP * a$GPP + b$GPP, 
                                    Reco = Reco * a$Reco + b$Reco, 
                                    Reco_sd = Reco * cv)

    # reform dataframe back to rasterLayer
    output.brk <- rasterFromXYZ(output.df[, c('x', 'y', 'Reco', 'Reco_sd')])
    crs(output.brk) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    names(output.brk) <- list('reco.mean', 'reco.sd')

    cat('predReco_loop_keras(): END of this function of predicting Reco at 500m')
    return(output.brk)
}   

# end of function

