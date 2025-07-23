# This code is used to calculate the percentage tree cover of each pixel over 
# the city of Toronto. The output is used in the shoreline correction
# for downscaled Solar Induced Fluorescence

# Parts of the code marked as *** are portions the user should change (e.g. paths)

memory.limit(size=5e5)

library(raster)
library(ncdf4)
library(ggplot2)

# *** CHANGE PATH for SIF data ***
sif_dir <-"C:/Users/kitty/Documents/Research/SIF/SMUrF/data/downscaled_CSIF/TROPOMI_CSIF_combined_med/V061/2018/V3/"
sif_files <- list.files(sif_dir)

# *** CHANGE PATH & FILE NAME for Toronto tree cover data ***
toronto_trees_high_res<-raster('C:/Users/kitty/Documents/ArcGIS/Projects/Toronto_canopy_cover/Toronto_tree_cover_high_res.tif')

sif_crop0 <- crop(raster(paste0(sif_dir,sif_files[20])),toronto_trees_high_res)

#-------------------------------------------------------------------------#
### *** Uncomment if first time running: ***

#toronto_trees_high_res[toronto_trees_high_res==0]<-NA
#toronto_trees_high_res[toronto_trees_high_res>1]<-0
#plot(toronto_trees_high_res)

#tree_resamp <- resample(toronto_trees_high_res,sif_crop0)
#plot(tree_resamp)
#tree_resamp[tree_resamp<0]<- -999
#tree_resamp[is.na(tree_resamp)] <- -999

#writeRaster(tree_resamp,paste0('C:/Users/kitty/Documents/ArcGIS/Projects/Toroonto_canopy_cover/Toronto_tree_cover_percent_rasmpled.nc'), overwrite=TRUE, varname="tree_cover", xname="lon", yname="lat")

## End of uncomment
#-----------------------------------------------------------------------#

# Load in tree cover percent (resampled to same resolution as downscaled SIF)
# *** Change Path ***
tree_resamp <- raster(paste0('C:/Users/kitty/Documents/ArcGIS/Projects/Toroonto_canopy_cover/Toronto_tree_cover_percent_rasmpled.nc'),varname="tree_cover")

tree_df <- as.data.frame(rasterToPoints(tree_resamp))
names(tree_df)[3]<-'tree_fraction'
tree_df$tree_fraction[tree_df$tree_fraction== -999] <- NA

#Load in land cover data
# *** Change Path ***
lc <- raster("C:/Users/kitty/Documents/Research/SIF/SMUrF/data/MCD12Q1/MCD12Q1.061_LC_Type1_doy2018001_aid0001.tif", varname = 'LC_Type1')
lc_crop <- crop(lc,tree_resamp) 
#plot(lc_crop) #*** Optional, uncomment for visualization

# Identify pixels over water
lc_df <- as.data.frame(rasterToPoints(lc_crop-sif_crop0*0))
names(lc_df)[3] <- "LandCover"
lc_df$LandCover[lc_df$LandCover<16] <- 0
lc_df$LandCover[lc_df$LandCover>=16] <- 1

# *** Uncomment Optional for visualization of water cover map: ***
#ggplot(lc_df,aes(x=x,y=y))+geom_tile(aes(fill=LandCover,width=1/240,
#        height=1/240))+ggtitle('Water, 2018')+labs(x=expression(longitude^0)
#          , y = expression(latitude^o))
# ---- End of optional ----

# remove negative SIF values & set sif to fill values when tree cover is fill
sif_crop0[sif_crop0<0] <- 0
sif_crop0[tree_resamp == -999]<- -999

sif_df0 <- as.data.frame(rasterToPoints(sif_crop0))
names(sif_df0)[3] <- "SIF"
sif_df0$SIF[sif_df0$SIF== -999] <- NA

sif_df0 <- cbind(sif_df0,lc_df$LandCover)
names(sif_df0)[4] <- "LandCover"

#set pixels over water to fill value
sif_df0$dist_to_wtr <- sif_df0$x*0-999

wtr_x <- sif_df0$x[sif_df0$LandCover==1]
wtr_y <- sif_df0$y[sif_df0$LandCover==1]

#Calculate the distance to the water
for(i in 1:length(sif_df0$SIF)){
  if(is.na(sif_df0$SIF[i])){
    sif_df0$dist_to_wtr[i] <- NA
  }else{
    xval <- sif_df0$x[i]
    yval <- sif_df0$y[i]
    sif_df0$dist_to_wtr[i] <- min(sqrt((wtr_x-xval)^2+(wtr_y-yval)^2))
  }
}

# *** Optional uncomment for visualisation of distance to water ***
#ggplot(sif_df0,aes(x=x,y=y))+geom_tile(aes(fill=dist_to_wtr<0.012,width=1/240,
#      height=1/240))+ggtitle('Distance to Water')+labs(x=expression(longitude^0)
#        , y = expression(latitude^o))
# --- End of optional ---


for (f in sif_files[1:46]){ # *** may need to subset files if other files in directory
  DoY <- substr(f,35,37)
  sif_dat<-raster(paste0(sif_dir,f))
  sif_crop <- crop(sif_dat,tree_resamp)
  #plot(sif_crop) #uncomment to plot original sif data
  sif_crop[sif_crop<0] <- 0
  sif_crop[tree_resamp == -999]<- -999
  
  sif_df <- as.data.frame(rasterToPoints(sif_crop))
  names(sif_df)[3] <- "SIF"
  sif_df$SIF[sif_df$SIF== -999] <- NA
  
  sif_tree_df<- cbind(sif_df,tree_df$tree_fraction)
  names(sif_tree_df)[4] <- 'tree_fraction'
  
  # Fit all Toronto sif data (including shoreline) to tree cover percentage
  fit <- lm(sif_tree_df$SIF ~ sif_tree_df$tree_fraction)
  m <- unname(coef(fit)[2])/100
  b <- unname(coef(fit)[1])
  R <- cor(sif_tree_df$SIF,sif_tree_df$tree_fraction,use='complete.obs')
  
  #Uncomment to visualize fit
  #print(ggplot(sif_tree_df,aes(x=tree_fraction*100,y=SIF))+geom_point()+labs(x=
  #        'Tree Fraction (%)', y='SIF (umol/m2/s)')+ggtitle(paste0(
  #          'SIF vs. tree cover, 500m, with shoreline, DoY = ',DoY,', 2018'))+geom_smooth(
  #            method=lm, se=FALSE, formula = y~x)+geom_text(x=18, y=max(sif_tree_df$SIF,na.rm=TRUE)+0.05,label=
  #              paste0(round(m,digits=5),' x + ',round(b,digits=5),
  #                ', R2 = ',round(R^2,digits=3)))+ylim(min(sif_tree_df$SIF,na.rm=TRUE)-0.01,max(sif_tree_df$SIF,na.rm=TRUE)+0.05))
  
  sif_tree_df$tree_ratio <- sif_tree_df$SIF/sif_tree_df$tree_fraction
  sif_tree_df$tree_ratio[sif_tree_df$tree_fraction==0]<-NA
  #sif_tree_df$ratio[sif_tree_df$ratio>500]<-NA
  
  #Uncomment to visualize ratio of sif to tree cover
  #ggplot(sif_tree_df,aes(x=x,y=y))+geom_tile(aes(fill=tree_ratio,width=1/240,
  #        height=1/240))+scale_fill_gradient2(limits=c(0,20))+ggtitle(paste0(
  #          'SIF/Tree Fraction, DoY = ',DoY,', 2018'))+labs(x=expression(longitude^0)
  #            , y = expression(latitude^o))
  
  #print(median(sif_tree_df$tree_ratio,na.rm=TRUE))
  #print(mean(sif_tree_df$tree_ratio,na.rm=TRUE))
  
  sif_tree_df$dist_to_wtr <- sif_df0$dist_to_wtr
  
  #Uncomment to visualize distance to water
  ##ggplot(sif_tree_df,aes(x=dist_to_wtr,y=tree_ratio))+geom_point()+labs(x=
  ##        'Distance to Water (o)', y='SIF/tree cover')+ggtitle(
  ##          'SIF-tree cover ratio vs. distance to lake')+ylim(0,5)
  
  df2 <- aggregate(sif_tree_df, by=list(cut(sif_tree_df$dist_to_wtr,seq(0,0.18,0.001))),median)
  
  # Uncomment to visualize sif-tree cover ratio vs. distance to water
  #print(ggplot(df2,aes(x=dist_to_wtr,y=tree_ratio))+geom_point()+labs(x=
  #        'Binned Distance to Water (o)', y='SIF/tree cover')+ggtitle(
  #        'Binned SIF-tree cover ratio vs. distance to lake')+
  #         ylim(0,max(df2$tree_ratio,na.rm=TRUE)+0.05)+geom_vline(
  #         xintercept=0.012)+geom_vline(xintercept = 0.034)+geom_text(
  #         x=0.004,y=max(df2$tree_ratio,na.rm=TRUE),label='x=0.012')+
  #         geom_text(x=0.0251,y=max(df2$tree_ratio,na.rm=TRUE),label='x=0.034'))
  
  sif_tree_df$SIF_no_wtr <- sif_tree_df$SIF
  sif_tree_df$SIF_no_wtr[sif_tree_df$dist_to_wtr<0.012] <- NA
  
  #remove oversampled data since it is on a different scale than the rest of the data
  num_uni <- sort(table(sif_tree_df$SIF_no_wtr),decreasing=TRUE)
  oversamp <- as.numeric(names(num_uni[num_uni>25]))
  
  sif_tree_df$fit_SIF <- sif_tree_df$SIF_no_wtr
  
  if(length(oversamp)>0){
    for(o in oversamp){
      sif_tree_df$fit_SIF[round(sif_tree_df$fit_SIF,digits=15)==round(o,digits=15) & sif_tree_df$fit_SIF!=0]<- NA
    }
  }
  
  ##ggplot(sif_tree_df,aes(x=x,y=y))+geom_tile(aes(fill=fit_SIF,width=1/240,
  ##        height=1/240))+ggtitle(paste0('SIF, DoY = ',DoY,', 2018'))+
  ##          labs(x=expression(longitude^0), y = expression(latitude^o))
  
  #sif_tree_df$SIF_no_wtr[sif_tree_df$SIF_no_wtr>0.24]<-NA
  
  # fit sif data to tree fraction with shoreline and oversampled sif removed
  fit <- lm(sif_tree_df$fit_SIF ~ sif_tree_df$tree_fraction)
  m <- unname(coef(fit)[2])/100
  b <- unname(coef(fit)[1])
  R <- cor(sif_tree_df$SIF_no_wtr,sif_tree_df$tree_fraction,use='complete.obs')
  
  #visualize the fit
  print(ggplot(sif_tree_df,aes(x=tree_fraction*100,y=fit_SIF))+geom_point()+labs(x=
        'Tree Fraction (%)', y='SIF (umol/m2/s)')+ggtitle(paste0(
            'SIF vs. tree cover, 500m, without shoreline and oversampling, DoY = '
              ,DoY,', 2018'))+geom_smooth(method=lm, se=FALSE, formula = y~x)+
                geom_text(x=18, y=max(sif_tree_df$SIF,na.rm=TRUE)+0.05,
                  label=paste0(round(m,digits=5),' x + ',round(b,digits=5),
                    ', R2 = ',round(R^2,digits=3)))+
                      ylim(y=min(sif_tree_df$SIF,na.rm=TRUE)-0.01,y=max(sif_tree_df$SIF,na.rm=TRUE)+0.05))
  
  pval <- t.test(x=sif_tree_df$tree_fraction, y=sif_tree_df$fit_SIF)[3]
  
  # If the correlation is non-negligible replace sif data near the water with tree-cover correction
  if(is.finite(R) & R^2>0.1 & pval$p.value<0.05){
    sif_tree_df$LandCover <- lc_df$LandCover
    
    sif_tree_df$SIF_corr <- sif_tree_df$SIF_no_wtr
    sif_tree_df$SIF_w <- (sif_tree_df$dist_to_wtr)/0.012
    sif_tree_df$SIF_w[!is.na(sif_tree_df$SIF_no_wtr)] <- 1
    sif_tree_df$tree_w <- (0.012-sif_tree_df$dist_to_wtr)/0.012
    sif_tree_df$tree_w[!is.na(sif_tree_df$SIF_no_wtr)] <- 0
    sif_tree_df$SIF_corr_w <- sif_tree_df$SIF_no_wtr
    sif_tree_df$SIF_corr_w[is.na(sif_tree_df$SIF_no_wtr)] <- (sif_tree_df$SIF[is.na(sif_tree_df$SIF_no_wtr)]*sif_tree_df$SIF_w[is.na(sif_tree_df$SIF_no_wtr)]+(sif_tree_df$tree_fraction[is.na(sif_tree_df$SIF_no_wtr)]*100*m +b)*sif_tree_df$tree_w[is.na(sif_tree_df$SIF_no_wtr)])/(sif_tree_df$SIF_w[is.na(sif_tree_df$SIF_no_wtr)]+sif_tree_df$tree_w[is.na(sif_tree_df$SIF_no_wtr)])
    sif_tree_df$SIF_corr_w[!is.na(sif_tree_df$tree_fraction) & sif_tree_df$tree_fraction<0.05] <- sif_tree_df$SIF[!is.na(sif_tree_df$tree_fraction) & sif_tree_df$tree_fraction<0.05]
    sif_tree_df$SIF_corr_w[sif_tree_df$LandCover==1] <- NA
    
    sif_tree_df$SIF_corr[is.na(sif_tree_df$SIF_no_wtr)] <- sif_tree_df$tree_fraction[is.na(sif_tree_df$SIF_no_wtr)]*100*m +b
    sif_tree_df$SIF_corr[!is.na(sif_tree_df$tree_fraction) & sif_tree_df$tree_fraction<0.05] <- sif_tree_df$SIF[!is.na(sif_tree_df$tree_fraction) & sif_tree_df$tree_fraction<0.05]
    sif_tree_df$SIF_corr[sif_tree_df$LandCover==1] <- NA
    
    # Uncomment to plot tree-cover corrected SIF vs. tree fraction
    #ggplot(sif_tree_df,aes(x=tree_fraction*100,y=SIF_corr))+geom_point()+labs(x=
    #        'Tree Fraction (%)', y='Corrected SIF (umol/m2/s)')+ggtitle(paste0(
    #          'Corrected SIF vs. tree cover, DoY = ',DoY,', 2018')) #+geom_smooth(
    #            #method=lm, se=FALSE, formula = y~x)#+geom_text(x=18, y=0.86,label=
    #              #paste0(round(m,digits=5),' x + ',round(b,digits=5),
    #                #', R2 = ',round(R^2,digits=3)))+ylim(-0.05,0.9)
    
    #*** Uncomment to plot the tree-cover corrected sif
    #print(ggplot(sif_tree_df,aes(x=x,y=y))+geom_tile(aes(fill=SIF_corr_w,
    #        width=1/240,height=1/240))+
    #        ggtitle(paste0('Shoreline & Oversample Weigthed Corrected SIF, DoY = ',
    #        DoY,', 2018'))+labs(x=expression(longitude^0),
    #        y = expression(latitude^o))+coord_equal())
    
    corr_sif_rstr <- rasterFromXYZ(cbind(sif_tree_df$x,sif_tree_df$y,sif_tree_df$SIF_corr_w))
    crs(corr_sif_rstr) <- crs(sif_crop0)
    
    #uncomment to visualize rasterized tree-cover corrected sif
    #plot(corr_sif_rstr) #same as plot above
    
    # Replace SIF data in Toronto with the tree cover-corrected SIF
    corr_sif_ext <- extend(corr_sif_rstr,sif_dat)
    corr_sif_ext[is.na(corr_sif_ext)] <- sif_dat[is.na(corr_sif_ext)]
    corr_df <- as.data.frame(rasterToPoints(corr_sif_ext))
    #Optional: plot the tree-corrected SIF over the entire domain
    #ggplot(corr_df,aes(x=x,y=y))+geom_tile(aes(fill=layer,width=1/240,height=1/240))+xlim(min(sif_tree_df$x),max(sif_tree_df$x))+ylim(min(sif_tree_df$y),max(sif_tree_df$y))+scale_fill_gradient2(limits=c(0,1))+coord_equal()
    
    writeRaster(corr_sif_ext,paste0(sif_dir,'downscaled_v061_TROPO_CSIF_shore_weighted_corrected_8d_2018',DoY,'.nc'), overwrite=TRUE, varname="daily_sif", varunit="umol m-2 s-1", xname="lon", yname="lat")
    print(paste0('DoY ',DoY,' corrected, R2 = ',round(R^2,4),', p-value = ',pval))
  }else{
    print(paste0('DoY ',DoY,' weak correlation no adjustment, R2 = ',round(R^2,4),', p-value = ',pval))
  }
}
