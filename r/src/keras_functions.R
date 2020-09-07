# all packges and subroutines needed for NN model in Keras

library(keras); library(reticulate); library(tensorflow); library(lmodel2)

# ---------------------------------------------------------------------------- #
# Display training progress by printing a single dot for each completed epoch.
print_dot_callback <- keras::callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 80 == 0) cat('\n')
    cat('.')
  }
)    


# ---------------------------------------------------------------------------- #
# 0. create model 
build_model <- function(n1 = 32, n2 = 8, act_func = c('relu', 'sigmoid')[1], ncol = 3) {
  
  model <- keras::keras_model_sequential() %>%
            layer_dense(units = n1, activation = act_func, input_shape = ncol) %>%
            layer_dense(units = n2, activation = act_func) %>%
            layer_dense(units = 1)
  
  model %>% compile(loss = 'mse', optimizer = optimizer_rmsprop(),
                    metrics = list('mean_absolute_error'))
  
  model
}   # end function build_model


# ---------------------------------------------------------------------------- #
prep_nn <- function(data.all, sel.igbp, f = c('TA_F_MDS', 'TS_F_MDS_1', 'GPP'), 
                    a, b, shuffleTF = F) {
    
    if (sel.igbp != 'ALL') sel.df <- data.all %>% filter(igbp == sel.igbp)
    else sel.df <- data.all 

    sel.df <- sel.df %>% select_if(is.numeric) %>% arrange(TIMESTAMP) %>% na.omit()

    # if indicated, shuffling before training
    if (shuffleTF) sel.df <- sel.df[sample(nrow(sel.df)), ]   # if shuffle

    # --------------------------------------------------------------------------- #
    # 1. normalize features, by min and max
    # select features
    sel.df$TA_norm <- (sel.df[, colnames(sel.df) == f[1]] - b$TA) / a$TA
    sel.df$TS_norm <- (sel.df[, colnames(sel.df) == f[2]] - b$TS) / a$TS
    sel.df$GPP_norm <- (sel.df[, colnames(sel.df) == f[3]] - b$GPP) / a$GPP
    sel.df$Reco_norm <- (sel.df$Reco - b$Reco) / a$Reco
    sel.df$row.indx <- seq(1, max(nrow(sel.df)))

    # % for training data, % for testing data
    bins <- trunc(seq(1, max(sel.df$row.indx), length = 6))
    holdouts <- sort(as.numeric(unlist(sample_n(data.frame(seq(1, 5)), 1))))
    sel.indx <- NULL; for (h in holdouts) sel.indx <- c(sel.indx, bins[h] : bins[(h + 1)])

    test.df  <- sel.df %>% filter(row.indx %in% sel.indx)
    train.df <- sel.df %>% filter(!row.indx %in% test.df$row.indx) %>% arrange(row.indx) 

    # get ready for the training and testing data and labels
    train_data <- as.matrix(train.df %>% dplyr::select('TA_norm', 'TS_norm', 'GPP_norm'))
    test_data <- as.matrix(test.df %>% dplyr::select('TA_norm', 'TS_norm', 'GPP_norm'))
    train_labels <- train.df$Reco_norm
    test_labels <- test.df$Reco_norm

    prep_list <- list(train_data = train_data, train_labels = train_labels, 
                      test_data = test_data, test_labels = test_labels)
    prep_list
}   # end function 



# ---------------------------------------------------------------------------- #
plot_nn_ds <- function(pred_y = pred_y, test_y = test_y, 
                       feature = c('FLUXNET', 'ERA5')[1], sel.igbp) {
    
    data_y <- data.frame(pred_y = pred_y, test_y = test_y)

    # calculate all errors in regular scales
    test_rmse <- sqrt(mean((pred_y - test_y)^2, na.rm = T))
    test_bias <- mean(pred_y - test_y, na.rm = T)
    test_cv   <- test_rmse / mean(test_y, na.rm = T)
    test_cor  <- cor(pred_y, test_y)
    test_lm  <- lmodel2(pred_y ~ test_y)$regression.results[1, ]; print(test_lm$Slope)
    test_str <- data.frame(s = c(1, test_lm$Slope), i = c(0, test_lm$Intercept), 
                          fac = c('1:1', as.character(test_lm$Method)))

    # ----
    ## plotting
    sz <- 0.6
    max.xy <- max(test_y, pred_y, na.rm = T)
    min.xy <- min(test_y, pred_y, na.rm = T)
    cat(paste0('max observed/predicted Reco: ', signif(max.xy, 3), 
            '; min observed/predicted Reco: ', signif(min.xy, 3), '\n'))

    # plot scatter points with density
    breaks <- seq(-6, 25, 2); font.size = rel(0.9)
    t1 <- ggplot() + theme_classic() + 
        scale_x_continuous(breaks = breaks, labels = breaks, 
                            limits = c(floor(min.xy), ceiling(max.xy))) + 
        scale_y_continuous(breaks = breaks, labels = breaks, 
                            limits = c(floor(min.xy), ceiling(max.xy))) + 
        geom_hex(data = data_y, aes(x = test_y, y = pred_y), bins = 50, size = sz) + 
        scale_fill_gradientn(colours = rev(terrain.colors(20))[-1], 
                             name = 'COUNT', trans = 'log10') + 
        geom_abline(data = test_str, aes(slope = s, intercept = i, linetype = fac), 
                    colour = 'gray30') + 
        scale_linetype_discrete(name = NULL) + 
        guides(linetype = guide_legend(nrow = 2)) +
        labs(x = 'Observed Reco', y = 'Predicted Reco') + 
        labs(title = paste('NN for Reco based on', feature, 'for', sel.igbp, 
                            '\nRMSE:', signif(test_rmse, 3),
                            '; Mean bias:', signif(test_bias, 2), 
                            '; CV:', signif(test_cv, 3), 
                            '; r:', signif(test_cor, 3))) + 
        theme(legend.position = c(0.13, 0.8), legend.key.width = unit(0.4, 'cm'),
                legend.text = element_text(size = font.size),
                legend.key = element_blank(), legend.key.height = unit(0.3, 'cm'),
                axis.title.y = element_text(size = font.size, angle = 90),
                axis.title.x = element_text(size = font.size, angle = 0),
                axis.text = element_text(size = font.size),
                axis.ticks = element_line(size = font.size),
                title = element_text(size = font.size),
                strip.text = element_text(size = font.size))
    t1
}
