#' a funcion inside predReco() to predict Reco based on pre-trained NN
#' @author DW, 08/01/2019

#' add biomes-specific Reco uncertainty, DW, 08/02/2019

#' update by DW:
#' 09/10/2019, use updated Land Cover tif generated from main_script_GPP.r
#' 09/19/2019, increase the Tair for urban by 2 deg C. 
#' 12/02/2019, for a few grid cells near water bodies whose temps are NA, 
#'             use a simple linear relation between GPP and Reco, ~ about half
#' 03/02/2020, instead of gap-filled land cover, reuse vegetated fractions. 

# ---------------------------------------------------------------------------- #
predReco.loopv2 <- function(prep.list, reg.name, timestr, skip.igbp, nn.dir, 
                            nn.pattern, nn.reg = c('conus', 'global')[2],
                            train.feature) {
    
    input.brk <- prep.list$input.brk 
    lc.rt     <- prep.list$lc.rt
    fveg.stk  <- prep.list$fveg.stk

    # predict Reco for each land cover type
    # treat CRO/NVM the same as CRO
    igbp.info <- get.igbp.col() %>% filter(!name %in% c(skip.igbp, 'CRO/NVM', 'WGT'))

    # before prediction, normalize all fields first 
    
    # search for all available NN models and get CV for Reco error
    nn.file <- list.files(nn.dir, nn.pattern, full.names = T)
    cv.df <- get.reco.cv(train.feature, reg = nn.reg)

    # no needs for 'WAT', 'BAR', 'SI'
    for (i in 1 : nrow(igbp.info)) {

        tmp.name <- igbp.info[i, 'name']
        tmp.cv <- cv.df %>% dplyr::filter(biomes == tmp.name)
        cat(paste0('predReco.loopv2(): working on predicting ', tmp.name, ';'))

        if (tmp.name == 'URB') {
            
            urban.gra.reco <- predReco.nnv2(input.brk, lc.rt, nn.file[grepl('_GRA', nn.file)], tmp.name)
            urban.dbf.reco <- predReco.nnv2(input.brk, lc.rt, nn.file[grepl('_DBF', nn.file)], tmp.name)
            urban.ebf.reco <- predReco.nnv2(input.brk, lc.rt, nn.file[grepl('_EBF', nn.file)], tmp.name)
            urban.enf.reco <- predReco.nnv2(input.brk, lc.rt, nn.file[grepl('_ENF', nn.file)], tmp.name)
            gc() 

            # get fractional uncertainty for GRA, DBF, EBF, and ENF 
            gra.cv <- cv.df[cv.df$biomes == 'GRA', 'cv']
            dbf.cv <- cv.df[cv.df$biomes == 'DBF', 'cv']
            ebf.cv <- cv.df[cv.df$biomes == 'EBF', 'cv']
            enf.cv <- cv.df[cv.df$biomes == 'ENF', 'cv']

            #' weight different Reco based on @param fveg.stk
            tmp.reco.mean <- urban.gra.reco * fveg.stk$GRA + urban.dbf.reco * fveg.stk$DBF + 
                             urban.ebf.reco * fveg.stk$EBF + urban.enf.reco * fveg.stk$ENF

            tmp.reco.sd <- sqrt(gra.cv ^ 2 * fveg.stk$GRA + dbf.cv ^ 2 * fveg.stk$DBF + 
                                ebf.cv ^ 2 * fveg.stk$EBF + enf.cv ^ 2 * fveg.stk$ENF) * tmp.reco.mean

        } else {
            tmp.nn.file   <- nn.file[grepl(paste0('_', tmp.name), nn.file)]
            tmp.reco.mean <- predReco.nnv2(input.brk, lc.rt, tmp.nn.file, tmp.name); gc()

            # also assign error CV to each grid 
            cat('assigning grid-level errors...\n')
            tmp.reco.sd <- tmp.reco.mean * tmp.cv$cv # CV in fractional uncertainty form

        }   # end if

        if (i == 1) { 
            reco.mean <- tmp.reco.mean; reco.sd <- tmp.reco.sd 
        } else {
            reco.mean <- raster::merge(reco.mean, tmp.reco.mean)
            reco.sd <- raster::merge(reco.sd, tmp.reco.sd)
        }   # end if

    }   # end for i 

    pred.stk <- stack(reco.mean, reco.sd)
    names(pred.stk) <- c('reco.mean', 'reco.sd')

    cat('predRecp.loopv2(): END of this function of predicting Reco at 500m')
    return(pred.stk)
}   
# end of predReco.loop()



# ---------------------------------------------------------------------------- #
#' inner function inside @function predReco.loopv2() starting at ~line 64
# ---------------------------------------------------------------------------- #
predReco.nnv2 <- function(input.brk, lc.rt, tmp.nn.file, tmp.name) {
    
    # read NN models and modify variable names
    nn.list <- readRDS(tmp.nn.file)
    tmp.nn  <- nn.list$nn

    # because trained model may have different names, 
    # e.g., Tair.era5, Tsoil.era5 and GPP.csif
    # force them to be Tair, Tsoil and GPP, so that we can do further prediction
    rev.nn <- tmp.nn; rev.var <- rev.nn$model.list$variables
    
    # if training model parameters == GPP.csif, Tair.* or Tsoil.*, 
    # force them to be GPP, Tair and Tsoil, in order to match the column names of df
    rev.var[grepl('GPP', rev.var)] <- 'GPP'
    rev.var[grepl('TA', rev.var)] <- 'TA'
    rev.var[grepl('TS', rev.var)] <- 'TS'
    rev.nn$model.list$variables <- rev.var
    dimnames(rev.nn$covariate)[[2]] <- rev.var

    # grab IGBP info and select land cover map for the desired land types
    igbp.info <- get.igbp.col()
    tmp.val   <- igbp.info[igbp.info$name == tmp.name, 'val']
    if (tmp.name == 'CRO') 
        tmp.val <- c(tmp.val, igbp.info[igbp.info$name == 'CRO/NVM', 'val'])

    tmp.lc.rt <- lc.rt
    tmp.lc.rt[!tmp.lc.rt %in% tmp.val] <- NA

    # select input variables for a certain land cover type
    x <- raster::mask(input.brk, tmp.lc.rt)

    if ('coef' %in% names(nn.list)) {

        cat('predReco.nnv2(): normalizing variables; ')
        a <- nn.list$coef$a; b <- nn.list$coef$b 
        x$TA  <- (x$TA - b[, grepl('TA', names(b))]) / a[, grepl('TA', names(a))]
        x$TS  <- (x$TS - b[, grepl('TS', names(b))]) / a[, grepl('TS', names(a))]
        x$GPP <- (x$GPP - b[, grepl('GPP', names(b))]) / a[, grepl('GPP', names(a))]

        cat('predicting variables...\n')
        x.Reco <- raster::predict(x, rev.nn)
        x.Reco <- x.Reco * a$Reco + b$Reco
        cat(paste('max Reco for', tmp.name, ':', signif(maxValue(x.Reco), 3), '\n'))

    } else x.Reco <- raster::predict(x, rev.nn)   # end if

    return(x.Reco)
}   # end function



# for sanity check 
if (F) {
    
    # easternUS cities
    e1 <- extent(-72, -69.5, 41, 43.5)
    e2 <- extent(-87, -85, 39, 41)
    e3 <- extent(-76, -74, 39, 41)
    e4 <- extent(-81.4, -80.4, 34.8, 35.8)
    e5 <- extent(-89, -87, 41.3, 43.3)
    e6 <- extent(-81, -79.5, 25, 27)
    sc <- c(1.4, 1.3, 1.3, 1.3, 1.2, 1.2)
    lab <- c('Boston', 'Indy', 'Philly-NY', 'Charlotte', 'Chicago-Milwaukee', 'Miami')
    #lab <- 'Cairo'; e1 <- extent(28, 33, 27, 32)

    # ---- for veg fractions
    all.stk <- stack(input.brk$GPP, reco.mean, raster::merge(reco.mean, urban.reco))
    names(all.stk) <- c('GPP', 'Reco_Before', 'Reco_After')

    plot.reco <- function(stk, ext) {
        col <- rasterTheme(region = c('gray90', rev(terrain.colors(12))[-1]))
        g1 <- grid.arrange(levelplot(crop(stk$GPP, ext), margin = F, 
                           maxpixel = 1E6, par.settings = col, at = seq(0, 24, 1), 
                           xlab = 'LONGITUDE', ylab = 'LATITUDE'))
        r1 <- grid.arrange(levelplot(crop(stk[[c(2, 3)]], ext), layout = c(1, 2), 
                           maxpixel = 1E6, par.settings = col, at = seq(0, 15, 1), 
                           xlab = 'LONGITUDE', ylab = 'LATITUDE'))
        gr <- ggarrange(g1, r1, nrow = 2, heights = c(1.2, 2))
        return(gr)
    }   # end function

    f1 <- plot.reco(all.stk, e1)
    f2 <- plot.reco(all.stk, e2)
    f3 <- plot.reco(all.stk, e3)
    f4 <- plot.reco(all.stk, e4)
    f5 <- plot.reco(all.stk, e5)
    f6 <- plot.reco(all.stk, e6)

    ff <- ggarrange(f1, f2, f3, f4, f5, f6, nrow = 1, ncol = 6, labels = lab) 
    ff <- annotate_figure(ff, top = 'Reprojected GPP & estimated ecosystem respiration over 500m MODIS-based URBAN areas')
    ggsave(ff, filename = '../Fig2_gpp_reco_east.png', width = 18, height = 9)


    # print(uni.igbp)
    #[1] "WET"     "ENF"     "EBF"     "MF"      "WSAV"    "GRA"     "DBF"
    #[8] "SAV"     "DNF"     "CRO"     "OSHR"    "ENF/DBF" "EBF/DBF" "CSHR"
    pred.df <- input.df %>% mutate(Reco = NA, Reco.sd = NA) %>% filter(!is.na(name))  
    
    tmp.igbp <- c('CRO', 'DBF', 'ENF', 'GRA', 'SAV', 
                  'WET', 'CSHR', 'MF', 'EBF', 'DNF', 
                  'WSAV')[11]
    tmp.df <- pred.df %>% filter(name == tmp.igbp) #%>% mutate(lab = 'initial')
    #add.df <- tmp.df %>% mutate(TA = TA - 20, TS = TS - 20, lab = 'fake')
    #tmp.df <- rbind(tmp.df, add.df)

    tmp.df$Reco <- predReco.nn(nn.file, tmp.igbp, tmp.df)
    tmp.df$Reco.sd = NULL
    tmp.df <- tmp.df %>% na.omit()
    print(range(tmp.df$Reco, na.rm = T))
    print(range(tmp.df$GPP, na.rm = T))

    #if (nrow(tmp.df) > 1E4) tmp.df <- sample_n(tmp.df, 5000)

    sz <- 0.5
    plot(tmp.df$GPP, tmp.df$Reco)

    a1 <- ggplot(data = tmp.df) + theme_bw() + geom_hline(yintercept = 0, linetype = 2) +
          labs(x = 'Modeled TA', y = 'Modeled Reco', title = tmp.igbp) + 
          geom_point(aes(TA, Reco, colour = GPP), size = sz) + 
          scale_colour_gradientn(colors = rev(heat.colors(5)), name = 'Modeled GPP')

    s1 <- ggplot(data = tmp.df) + theme_bw() + geom_hline(yintercept = 0, linetype = 2) +
          labs(x = 'Modeled TS', y = 'Modeled Reco', title = tmp.igbp) +
          geom_point(aes(TS, Reco, colour = GPP), size = sz) + 
          scale_colour_gradientn(colors = rev(heat.colors(5)), name = 'Modeled GPP')

    as <- ggarrange(a1, s1, nrow = 2, common.legend = T)


}