library(RColorBrewer)
library(neuralnet)

# ------------------------------------------------------------------------ #
#' inner function inside @function predReco_loop_neuralnet()
# ------------------------------------------------------------------------ #
neuralnet.reco <- function(input.mk, tmp.nn.fn) {
    
    # read NN models and modify variable names
    if (!file.exists(tmp.nn.fn)) 
        stop(paste('neuralnet.reco(): NN model file', basename(tmp.nn.fn), 
                   'not found, please check...\n'))

    nn.list <- readRDS(tmp.nn.fn)
    tmp.nn  <- nn.list$nn

    # because trained model may have different names, e.g., TA_era5, TS_era5 and GPP_mean
    # force them to be TA, TS and GPP, so that we can do further prediction
    rev.nn <- tmp.nn; rev.var <- rev.nn$model.list$variables
    
    # if training model parameters == GPP.csif, Tair.* or Tsoil.*, 
    # force them to be GPP, Tair and Tsoil, in order to match the column names of df
    rev.var[grepl('GPP', rev.var)] <- 'GPP'
    rev.var[grepl('TA', rev.var)] <- 'TA'
    rev.var[grepl('TS', rev.var)] <- 'TS'
    rev.nn$model.list$variables <- rev.var
    dimnames(rev.nn$covariate)[[2]] <- rev.var

    # predict reco
    y.reco <- raster::predict(input.mk, rev.nn)

    # remove Reco outlier at 500 m 
    y.reco[y.reco < -0.1] <- 0
    y.reco[y.reco > 1.1] <- 1.0 

    return(y.reco)
}   # end function


# ------------------------------------------------------------------------ #
#' inner function inside @function predReco_loop_neuralnet(), see line 111 -
# ------------------------------------------------------------------------ #
neuralnet.reco.urban <- function(x, igbp.info, fveg.stk, nn.fns, nn.pattern) {

    y.gra <- neuralnet.reco(x, nn.fns[grepl('_GRA', nn.fns)])
    y.enf <- neuralnet.reco(x, nn.fns[grepl('_ENF', nn.fns)])
    y.mf  <- neuralnet.reco(x, nn.fns[grepl('_MF', nn.fns)])
    y.dbf <- neuralnet.reco(x, nn.fns[grepl('_DBF', nn.fns)])

    # get fractional uncertainty for GRA, DBF, EBF, and ENF 
    gra.cv <- igbp.info[igbp.info$name == 'GRA', 'cv']
    dbf.cv <- igbp.info[igbp.info$name == 'DBF', 'cv']
    enf.cv <- igbp.info[igbp.info$name == 'ENF', 'cv']
    mf.cv  <- igbp.info[igbp.info$name == 'MF', 'cv']

    if (nn.pattern == 'daymet_nldas') {
        y.shr  <- neuralnet.reco(x, nn.fns[grepl('_OSHR', nn.fns)])
        shr.cv <- igbp.info[igbp.info$name == 'OSHR', 'cv']

        # modify fveg.stk since no DNF and EBF is associated with CONUS 
        # when using daymet + NLDAS
        cat('neuralnet.reco.urban(): modifying vegetated fractions for urban when using Daymet + NLDAS...\n')
        fveg.sub <- subset(fveg.stk, which(!names(fveg.stk) %in% c('DNF', 'EBF')))
        fveg.mod <- fveg.sub / sum(fveg.sub); names(fveg.mod) <- names(fveg.sub)

        # weight Reco based on fveg.mod being proccessed during GPP estimates
        cat('neuralnet.reco.urban(): calculating weighted mean Reco for urban...\n')
        y <- y.gra * fveg.stk$GRA + y.shr * fveg.stk$OSHR + y.mf * fveg.stk$MF + 
             y.dbf * fveg.stk$DBF + y.enf * fveg.stk$ENF

        y.cv <- sqrt((gra.cv * fveg.stk$GRA) ^ 2 + (shr.cv * fveg.stk$OSHR) ^ 2 + 
                     (dbf.cv * fveg.stk$DBF) ^ 2 + (enf.cv * fveg.stk$ENF) ^ 2 + 
                     (mf.cv * fveg.stk$MF) ^ 2) 

    } else if (nn.pattern == 'era5') {    # if using ERA5

        y.shr <- neuralnet.reco(x, nn.fns[grepl('_SHR', nn.fns)])
        y.dnf <- neuralnet.reco(x, nn.fns[grepl('_DNF', nn.fns)])
        y.ebf <- neuralnet.reco(x, nn.fns[grepl('_EBF', nn.fns)])

        shr.cv <- igbp.info[igbp.info$name == 'SHR', 'cv']
        dnf.cv <- igbp.info[igbp.info$name == 'DNF', 'cv']
        ebf.cv <- igbp.info[igbp.info$name == 'EBF', 'cv']

        # weight Reco based on fveg.stk being proccessed during GPP estimates
        cat('nerualnet.reco.urban(): calculating weighted mean Reco for urban...\n')
        y <- y.gra * fveg.stk$GRA + y.shr * fveg.stk$OSHR + 
             y.dbf * fveg.stk$DBF + y.dnf * fveg.stk$DNF + 
             y.ebf * fveg.stk$EBF + y.enf * fveg.stk$ENF + y.mf * fveg.stk$MF

        y.cv <- sqrt((gra.cv * fveg.stk$GRA) ^ 2 + (shr.cv * fveg.stk$OSHR) ^ 2 + 
                     (dbf.cv * fveg.stk$DBF) ^ 2 + (dnf.cv * fveg.stk$DNF) ^ 2 + 
                     (ebf.cv * fveg.stk$EBF) ^ 2 + (enf.cv * fveg.stk$ENF) ^ 2 + 
                     (mf.cv * fveg.stk$MF) ^ 2) 
    } else {
        stop('neuralnet.reco.urban(): invalid nn.pattern...\n')
    
    }  # end if

    y.out <- stack(y, y.cv); names(y.out) <- c('y', 'y.cv')
    y.out
}   # end of nerualnet.reco.urban()



# --------------------------------------------------------------------------- #
# script to train Reco using NN, Dien Wu, 05/16/2019 
# add different features for predicting from training, DW, 08/01/2019
# always perform normalization before training, DW, 12/29/2019 
# --------------------------------------------------------------------------- #
train.reco <- function(data.all, sel.igbp, td = 0.8, hidden = c(32, 8), rep = 3, 
                       stepmax = 1000, overwriteTF = F, 
                       sample = c('shuffle', 'order')[1],   # 'order' for Chronological
                       reg = c('conus', 'global')[2], today = '20200401') {

    cat(paste('\n\n# --------- Working on', sel.igbp, '------------ #\n'))
    df <- data.all %>% filter(igbp == sel.igbp)
    if (sel.igbp == 'ALL') df <- data.all
    if ('RECO' %in% colnames(df)) df <- df %>% rename(Reco = RECO)

    uni.site <- unique(df$site); print(uni.site)
    num.site <- length(uni.site); print(num.site)
    note <- paste('for', num.site, sel.igbp, toupper(reg), 'sites'); print(nrow(df))

    # output resultant NN model into rds files 
    rds1 <- paste0('NN_model_Reco_tair_tsoil_gpp_', sel.igbp, '_fluxnet.rds')
    rds2 <- paste0('NN_model_Reco_tair_tsoil_sifgpp_', sel.igbp, '_era5.rds')
    rds3 <- paste0('NN_model_Reco_tair_tsoil_sifgpp_', sel.igbp, '_daymet+nldas.rds')

    #rds.file0 <- file.path(nn.path, 'train_daymet_nldas', rds0)
    rds.file1 <- file.path(nn.path, 'train_fluxnet', reg, rds1)
    rds.file2 <- file.path(nn.path, 'train_era5', reg, rds2)
    rds.file3 <- file.path(nn.path, 'train_daymet_nldas', reg, rds3)

    # separate training vs. testing data
    # remove outliers, very negative GPP, or low Tsoil but high Reco
    df <- df %>% filter(GPP_mean > min(df$GPP), TS_F_MDS_1 > 1 | Reco < 3) %>% 
                 arrange(TIMESTAMP) 
    
    # --------------------------------------------------------------------- #
    # plot raw relation first 
    # training features, should be found in colnames of df
    f1 <- c('TA_F_MDS', 'TS_F_MDS_1', 'GPP')     # all observed T + GPP
    f2 <- c('TA_era5', 'TS_era5', 'GPP_mean')    # both modeled T and GPP
    
    all.f <- c(f1, f2); nrow = 2; height = 8
    if (reg == 'conus') { 
        f3 <- c('TA_daymet', 'TS_nldas', 'GPP_mean')         # both modeled T and GPP
        all.f <- c(f1, f2, f3); nrow = 3; height = 12 
    }

    pp <- plot.reco.feature.tower(obs.df = df, features = all.f, nrow = nrow, ncol = length(f1))
    title <- paste0('Daily mean Reco vs. normalized X (2010-2014) ', note)
    pp <- annotate_figure(pp, top = text_grob(title))
    #ggsave(pp, filename = '../paper3/test_cumsum.png', width = 20, height = 10)

    # --------------------------------------------------------------------- #
    # reco.nn() returns a ggplot and automatically stores the NN model in rds format
    # with default rprop+ optimizer
    ### call reco.nn()
    nn11 <- reco.nn(df, f1, hidden, td, rep, stepmax, sample, note, overwriteTF, 
                    rds.file = rds.file1)  # FLUXNET GPP + T
    
    nn22 <- reco.nn(df, f2, hidden, td, rep, stepmax, sample, note, overwriteTF, 
                    rds.file = rds.file2)  # ERA5 temps

    if (reg == 'conus') 
        nn33 <- reco.nn(df, f3, hidden, td, rep, stepmax, sample, note, 
                        overwriteTF, rds.file = rds.file3)  # ERA5 temps

    # --------------------------------------------------------------------- #
    #nnn <- ggarrange(nn11, nn22, nrow = 2, common.legend = T, labels = c('d)', 'h)'))
    #if (reg == 'conus') nnn <- ggarrange(nn11, nn22, nn33, nrow = 3, common.legend = T)
    #pn <- ggarrange(pp, nnn, ncol = 2, widths = c(2.9, 1))
    #fn <- paste0('pred_train_reco_', sel.igbp, '_', reg, '_fluxnet_', today, '.png')
    #ggsave(pn, filename = file.path(nn.path, fn), width = 14, height = height)

    nn.list <- list(nn11 = nn11 %>% mutate(IGBP = sel.igbp), 
                    nn22 = nn22 %>% mutate(IGBP = sel.igbp))
    if (reg == 'conus') nn.list <- list(nn11 = nn11, nn22 = nn22, nn33 = nn33)

    return(nn.list)

}   # end of function



# --------------------------------------------------------------------------- #
# try ANN for training data
# --------------------------------------------------------------------------- # 
reco.nn <- function(df, train.feature = c('TA_F_MDS', 'TS_F_MDS_1', 'GPP'), 
                    hidden = c(8, 4), td = 2.5, rep = 3, stepmax = 8000, 
                    sample = c('shuffle', 'order')[2], note = NULL, 
                    overwriteTF = T, rds.file = NULL, font.size = rel(0.7)) {
    
    library(neuralnet); library(lmodel2); library(ggplot2)

    n <- names(df)
    f.train <- paste(n[n %in% train.feature], collapse = ' + ')
    f <- as.formula(paste('Reco ~', f.train))
    sel.df <- df %>% select_if(is.numeric)

    # perform MAX-MIN NORMALIZATION
    min.df <- data.frame('Reco' = -1, 'TA' = -35, 'TS' = -30, 'GPP' = -1)
    max.df <- data.frame('Reco' = 20, 'TA' = 40, 'TS' = 40, 'GPP' = 30)

    # calculate the coeffient to convert normalized scales back to regular scales
    # e.g., x = x_norm * a + b where x_norm fall within [0, 1]
    a <- max.df - min.df; b <- min.df 

    #sel.df <- as.data.frame(lapply(sel.df, normalize)) 
    #sel.df$TIMESTAMP <- df$TIMESTAMP
    sel.df <- sel.df %>% select_if(colnames(sel.df) %in% c('Reco', train.feature, 'TIMESTAMP'))
    
    # x_norm = (x - b) / a
    sel.df$TA_norm <- (sel.df[, colnames(sel.df) == train.feature[1]] - b$TA) / a$TA
    sel.df$TS_norm <- (sel.df[, colnames(sel.df) == train.feature[2]] - b$TS) / a$TS
    sel.df$GPP_norm <- (sel.df[, colnames(sel.df) == train.feature[3]] - b$GPP) / a$GPP
    sel.df$Reco_norm <- (sel.df$Reco - b$Reco) / a$Reco

    f.train <- paste('TA_norm + TS_norm + GPP_norm')
    f <- as.formula(paste('Reco_norm ~', f.train))

    cat('\n\n# --------------------------------------------- #\n')
    # order by time (not by sites) before training 
    sel.df <- sel.df %>% arrange(TIMESTAMP)

    # if indicated, shuffling before training
    if (sample == 'shuffle') sel.df <- sel.df[sample(nrow(sel.df)), ]
    sel.df$row.indx  <- seq(1, max(nrow(sel.df)))

    # 80% for training data, 20% for testing data
    bins <- trunc(seq(1, max(sel.df$row.indx), length = 6))
    k.holdout <- sample(5)[1]   # randomly pick a holdout
    sel.indx  <- seq(bins[k.holdout], bins[k.holdout + 1])
    test.df  <- sel.df %>% filter(row.indx %in% sel.indx)
    train.df <- sel.df %>% filter(!row.indx %in% test.df$row.indx) %>% arrange(row.indx) 
    
    #train.df <- sample_n(sel.df, nrow(sel.df) * 2 / 3) %>% arrange(row.indx) 
    #test.df  <- sel.df %>% filter(!sel.df$row.indx %in% train.df$row.indx) 

    if (file.exists(rds.file) & overwriteTF == F) { # found a rds file

        nn.list <- readRDS(rds.file)
        nn  <- nn.list$nn 
        td  <- nn.list$params$threshold
        rep <- nn.list$params$rep 
        stepmax <- nn.list$params$stepmax

        if ('coef' %in% names(nn.list)) {   # if normalized, read coef
            coef <- nn.list$coef 
            a <- coef$a 
            b <- coef$b
        }   # end if

    } else {    # no rds file found or need overwriting

        # fit neural network for ecosystem resp estimates
        cat('reco.nn(): Start training data...\n')
        # for debugging
        # 
        nn <- neuralnet(f, data = train.df, err.fct = 'sse', act.fct = 'tanh', 
                           hidden = hidden, rep = rep, threshold = td,
                           linear.output = T, stepmax = stepmax, 
                           lifesign = 'full')
                           
        cat('reco.nn(): Finish training data...\n')
    }  # end if file.exists


    # --------------------------------------------------------------------- #
    # predict data
    test.df$Reco_pred_norm  <- as.numeric(predict(nn, test.df))
    train.df$Reco_pred_norm <- as.numeric(predict(nn, train.df))

    # convert Reco and Reco_pred to their initial scales and calculate error stat
    # firstly rename all colnames 
    plot(test.df$Reco_norm, test.df$Reco_pred_norm)
    test.df <- test.df %>% mutate(Reco = Reco_norm * a$Reco + b$Reco, 
                                  Reco_pred = Reco_pred_norm * a$Reco + b$Reco)
    train.df <- train.df %>% mutate(Reco = Reco_norm * a$Reco + b$Reco, 
                                    Reco_pred = Reco_pred_norm * a$Reco + b$Reco)

    # --------------------------------------------------------------------- #
    rmse.test <- sqrt(mean((test.df$Reco_pred - test.df$Reco)^2, na.rm = T))
    bias.test <- mean(test.df$Reco_pred - test.df$Reco, na.rm = T)
    cv.test   <- rmse.test / mean(test.df$Reco, na.rm = T)
    cor.test  <- cor(test.df$Reco_pred, test.df$Reco)

    rmse.train <- sqrt(mean((train.df$Reco_pred - train.df$Reco)^2, na.rm = T))
    bias.train <- mean(train.df$Reco_pred - train.df$Reco, na.rm = T)
    cv.train   <- rmse.train / mean(train.df$Reco, na.rm = T)
    cor.train  <- cor(train.df$Reco_pred, train.df$Reco)

    cat(paste('reco.nn(): RMSE:', signif(rmse.test, 3), 
              '; Mean bias:', signif(bias.test, 3), 
              '; CV:', signif(cv.test, 3),  '; r:', signif(cor.test, 3), 
              'for testing data\n'))
    
    cat(paste('reco.nn(): RMSE:', signif(rmse.train, 3), 
              '; Mean bias:', signif(bias.train, 3), 
              '; CV:', signif(cv.train, 3),  '; r:', signif(cor.train, 3), 
              'for training data\n'))

    # store model
    if (!file.exists(rds.file) | overwriteTF) {  # store hyper parameters
        nn.list <- list(nn = nn, 
                        coef = list(a = a, b = b), 
                        params = list(hidden = hidden, rep = rep, threshold = td, 
                                      stepmax = stepmax, err.fct = 'sse', 
                                      act.fct = 'tanh', feature = train.feature), 
                        errors = list(rmse = rmse.test, bias = bias.test, cv = cv.test))
        saveRDS(nn.list, file = rds.file)
        cat(paste('reco.nn(): NN model stored in', rds.file, '\n\n'))
    }   # end if 

    # --------------------------------------------------------------------- #
    ## plotting
    sz <- 0.6
    max.xy <- max(test.df$Reco, test.df$Reco_pred, 
                  train.df$Reco, train.df$Reco_pred, na.rm = T)
    min.xy <- min(test.df$Reco, test.df$Reco_pred, 
                  train.df$Reco, train.df$Reco_pred, na.rm = T)
    cat(paste0('max observed/predicted Reco: ', signif(max.xy, 3), 
               '; min observed/predicted Reco: ', signif(min.xy, 3), '\n'))

    # fit a linear regression
    lm.fit.test <- lmodel2(Reco_pred ~ Reco, test.df)$regression.results[1, ]
    lm.fit.train <- lmodel2(Reco_pred ~ Reco, train.df)$regression.results[1, ]
    print(lm.fit.test$Slope)
    lm.test.df <- data.frame(s = c(1, lm.fit.test$Slope), 
                             i = c(0, lm.fit.test$Intercept), 
                             fac = c('1:1', as.character(lm.fit.test$Method)))
    lm.train.df <- data.frame(s = c(1, lm.fit.train$Slope), 
                              i = c(0, lm.fit.train$Intercept), 
                              fac = c('1:1', as.character(lm.fit.train$Method)))

    # plot scatter points with density
    breaks <- seq(-6, 25, 2)        
    t1 <- ggplot() + theme_classic() + 
          scale_x_continuous(breaks = breaks, labels = breaks, 
                             limits = c(-1, ceiling(max.xy))) + 
          scale_y_continuous(breaks = breaks, labels = breaks, 
                             limits = c(-1, ceiling(max.xy))) + 
          geom_hex(data = test.df, aes(x = Reco, y = Reco_pred), bins = 50, size = sz) + 
          scale_fill_gradientn(colours = c('gray92', brewer.pal(9, 'YlGnBu')),    # rev(terrain.colors(20))[-1]
                               name = 'COUNT', trans = 'log10') + 
          geom_abline(data = lm.test.df, aes(slope = s, intercept = i, linetype = fac), 
                      colour = 'gray30') + 
          scale_linetype_discrete(name = NULL) + guides(linetype = guide_legend(nrow = 2)) +
          labs(x = 'Observed Reco', y = 'Predicted Reco') + 
          labs(title = paste('NN for Reco ~', f.train, 
                             #'\nwith threshold of', td, '; rep of', rep, 
                             '\nRMSE:', signif(rmse.test, 3),
                             '; Mean bias:', signif(bias.test, 2), 
                             '; CV:', signif(cv.test, 3), 
                             '; r:', signif(cor.test, 3), '\nfor', note)) + 
          theme(legend.position = c(0.13, 0.8), legend.key.width = unit(0.4, 'cm'),
                legend.text = element_text(size = font.size),
                legend.key = element_blank(), legend.key.height = unit(0.3, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))

    #return(t1)
    return(test.df)
}



# --------------------------------------------------------------------------- #
# script to carry out K-fold X-validation in R, DW, 03/23/2020
# modified from reco.nn() 
# --------------------------------------------------------------------------- #
kfold.xval <- function(all.df, sel.igbp, k = 5, 
                       train.feature = c('TA_F_MDS', 'TS_F_MDS_1', 'GPP'), 
                       sample = c('shuffle', 'order')[2], rep = 2, threshold = 1) {
    
    library(neuralnet); library(lmodel2); library(ggplot2)

    # tested layer combinations
    hidden.vec <- list(t1 = c(64, 8), t2 = c(32, 8), t3 = c(16, 8), t4 = c(8, 8), 
                       t5 = c(64, 4), t6 = c(32, 4), t7 = c(16, 4), #t8 = c(8, 4), 
                       t10 = c(64, 16), t11 = c(32, 16), t12 = c(16, 16), 
                       #t13 = 64, t14 = 32, t15 = 16, t16 = 8, 
                       t17 = c(64, 32, 16), t18 = c(64, 16, 8), t19 = c(32, 16, 8), 
                       t20 = c(16, 8, 4)) 
    
    # select data
    df <- all.df %>% filter(igbp == sel.igbp); n <- names(df)
    f.train <- paste(n[n %in% train.feature], collapse = ' + ')
    f <- as.formula(paste('Reco ~', f.train))

    # normalization 
    sel.df <- df %>% select_if(is.numeric)

    # 1. select columns and perform MAX-MIN NORMALIZATION
    min.df <- data.frame('Reco' = 0, 'TA' = -35, 'TS' = -30, 'GPP' = -1)
    max.df <- data.frame('Reco' = 20, 'TA' = 40, 'TS' = 40, 'GPP' = 30)

    # calculate the coeffient to convert normalized scales back to regular scales
    # e.g., x = x_norm * a + b where x_norm fall within [0, 1]
    a <- max.df - min.df; b <- min.df  
    sel.df <- sel.df %>% select_if(colnames(sel.df) %in% c('Reco', train.feature, 'TIMESTAMP'))
    
    sel.df$TA_norm <- (sel.df[, colnames(sel.df) == train.feature[1]] + 30) / 70
    sel.df$TS_norm <- (sel.df[, colnames(sel.df) == train.feature[2]] + 30) / 70
    sel.df$GPP_norm <- (sel.df[, colnames(sel.df) == train.feature[3]] + 1) / 31
    sel.df$Reco_norm <- (sel.df$Reco - 0) / 20
    f.train <- paste('TA_norm + TS_norm + GPP_norm')
    f <- as.formula(paste('Reco_norm ~', f.train))

    # ------------------------------------------------------------------------ #
    cat('\n\n# -------------- Start K-fold cross validation ------------- #\n')
    sel.df$row.indx <- as.numeric(rownames(sel.df))

    # fit neural network for ecosystem resp estimates
    stat.all <- NULL
    for (h in 1: length(hidden.vec)) {

        cat('kfold.xval(): Start K-fold validation for different hidden layers...\n')
        hidden <- hidden.vec[[h]]
        
        # order by time (not by sites) before training 
        sel.df <- sel.df %>% arrange(TIMESTAMP)

        # if indicated, shuffling before training
        if (sample == 'shuffle') sel.df <- sel.df[sample(nrow(sel.df)), ]
        rownames(sel.df) <- seq(1, max(sel.df$row.indx))
        sel.df$row.indx  <- seq(1, max(sel.df$row.indx))

        for (k.holdout in 1 : k) {
            stat.df <- kfold(sel.df, f, k, hidden, k.holdout, a, b, rep, threshold) %>% 
                       mutate(igbp = sel.igbp)
                        
            stat.all <- rbind(stat.all, stat.df)
        }   # end for k.holdout

    }   # end for h 

    stat.all <- stat.all %>% 
                mutate(hidden = ifelse(is.na(hidden3), 
                                        paste(hidden1, hidden2), 
                                        paste(hidden1, hidden2, hidden3)))
    return(stat.all)

}   # end of function



# --------------------------------------------------------------------------- #
# kfold X-val inner functions
# --------------------------------------------------------------------------- #
kfold <- function(sel.df, f, k, hidden, k.holdout, a, b, rep = 1, threshold = 1) {
    
    # perform K-fold cross validation, K for # of holdouts
    cat(paste('kfold(): #', k.holdout, 'holdout\n'))

    bins <- trunc(seq(1, max(sel.df$row.indx), length = k + 1))
    sel.indx <- seq(bins[k.holdout], bins[k.holdout + 1])
    val.df <- sel.df %>% filter(row.indx %in% sel.indx)
    train.df <- sel.df %>% filter(!row.indx %in% val.df$row.indx) %>% arrange(row.indx) 
    
    # NN training
    stepmax = 2000#; if (length(hidden) == 1) stepmax = 2000
    nn <- neuralnet(f, data = train.df, err.fct = 'sse', act.fct = 'tanh', 
                       hidden = hidden, rep = rep, threshold = threshold, 
                       linear.output = T, stepmax = stepmax, lifesign = 'full')
    
    # predict data and convert Reco and Reco_pred to their initial scales 
    # and calculate error stat
    val.df <- val.df %>% mutate(Reco_pred_norm = as.numeric(predict(nn, val.df)),
                                Reco = Reco_norm * a$Reco + b$Reco, 
                                Reco_pred = Reco_pred_norm * a$Reco + b$Reco)

    train.df <- train.df %>% mutate(Reco_pred_norm = as.numeric(predict(nn, train.df)), 
                                    Reco = Reco_norm * a$Reco + b$Reco, 
                                    Reco_pred = Reco_pred_norm * a$Reco + b$Reco)
    plot(val.df$Reco_norm, val.df$Reco_pred_norm)

    # --------------------------------------------------------------------- #
    # for validation data set, 1 / k
    rmse.val <- sqrt(mean((val.df$Reco_pred - val.df$Reco)^2, na.rm = T))
    bias.val <- mean(val.df$Reco_pred - val.df$Reco, na.rm = T)
    cv.val   <- rmse.val / mean(val.df$Reco, na.rm = T)
    cor.val  <- cor(val.df$Reco_pred, val.df$Reco)

    # for training data set, 1 - 1 / k 
    rmse.train <- sqrt(mean((train.df$Reco_pred - train.df$Reco)^2, na.rm = T))
    bias.train <- mean(train.df$Reco_pred - train.df$Reco, na.rm = T)
    cv.train   <- rmse.train / mean(train.df$Reco, na.rm = T)
    cor.train  <- cor(train.df$Reco_pred, train.df$Reco)

    # save all results 
    stat.df <- data.frame(k.holdout = k.holdout, hidden1 = hidden[1], 
                          hidden2 = hidden[2], hidden3 = hidden[3], 
                          rmse.val = rmse.val, bias.val = bias.val, 
                          cv.val = cv.val, cor.val = cor.val, 
                          rmse.train = rmse.train, bias.train = bias.train, 
                          cv.train = cv.train, cor.train = cor.train)

    return(stat.df)
}
# end of function
