#' @function This script takes the given year and region and converts MODIS 
#' reflectance to EVI for use in adjusting respiration in urban areas 
#' 
#' @author Sabrina Madsen, 02/19/2024
#' 

mod_EVI <- function(mod_dir='C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/dataverse_files/MODIS_reflectance/MODIS_V061_GTA_AppEEARS_2021',
                    y=2021,
                    minlon = -125, # in deg, negative for west hemisphere
                    maxlon = -95,  # in deg
                    minlat = 30, 
                    maxlat = 50,
                    reg.name,
                    lc.path,
                    lc.pattern,
                    lc.max.yr,
                    evi.path,
                    evi.pattern
                    ){
  setwd(mod_dir)
  mod_files <- list.files()
  b01_files <- mod_files[grep('b01',mod_files)]
  b02_files <- mod_files[grep('b02',mod_files)]
  b03_files <- mod_files[grep('b03',mod_files)]
  qc_files <- mod_files[grep('qc_500m_doy',mod_files)]
  
  #import land cover
  reg.ext <- raster::extent(minlon, maxlon, minlat, maxlat) # regional extent
  lc.rt <- prep.mcd12(lc.path, lc.pattern, y, lc.max.yr, reg.name, reg.ext)
  
  ##create a raster
  EVI_raster_V061 <-raster()
  ##set the number of columns, rows, and extent
  EVI_raster_V061 <- brick(ncol=ncol(lc.rt), nrow=nrow(lc.rt), xmn=minlon, xmx=maxlon, ymn=minlat, ymx=maxlat,nl=380)
  res(EVI_raster_V061)
  ##check the number of cells is 344448 for the GTA, 383040 for Montreal/Ottawa
  ncell(EVI_raster_V061)
  names(EVI_raster_V061)<-paste('DOY',-3:376,sep="") #name the layers
  
  
  
  for (i in 1:length(b01_files)){
    
    b01_file <- stack(b01_files[i])
    b02_file <- stack(b02_files[i])
    b03_file <- stack(b03_files[i])

    qc_file <- stack(qc_files[i])
    #Reproject the reflectance data to the land cover data
    file1<-raster::projectRaster(b01_file, lc.rt)
    file2<-raster::projectRaster(b02_file, lc.rt)
    file3<-raster::projectRaster(b03_file, lc.rt)
    qc_file <- raster::projectRaster(qc_file,lc.rt)

    # APPLY SCALE FACTOR FOR GTA 2018 BUT NOT FOR OTHER YEARS (THE DATA ON 
    # APPEEARS CHANGED)
    # *** THIS WILL NEED TO BE CHANGED DEPENDING ON WHEN YOU DOWNLOADED THE 
    #     MODIS DATA !!!!! ***
    ## apply scale factors for reflectance data in the bands needed for EVI/LSWI calculation
    if (y==2018 & length(grep('Toronto',lc.rt))>0){
      file1 <- file1 * 0.0001 #red
      file2 <- file2 * 0.0001 #NIR
      file3 <- file3 * 0.0001 #blue
    }
    
    ok_qc=c(1073741824, 1073954817, 1075838976, 1075838977, 1075838979, 1076051969, 1075576832,
            1076625411, 1077149697, 1077362689, 1077411843, 1119879171, 1121976323,
            1123287043, 1123549187, 1128267777, 1130364929, 1130364931, 1131151363,
            1131675649, 1131937795, 1132462083, 1134559235, 1135345667, 1135869955,
            1136132099, 1811939331, 1814036483, 1815347203, 1815609347, 1858076675,
            1866465283, 1868562435, 1869873155, 1870135299, 1870659587, 1873543171,
            1874329603, 1946157057, 1946370049, 1948254209, 1948254211, 1949040643,
            1949564929, 1949827075, 1992294403, 1995964419, 2000683009, 2000896001,
            2002780161, 2002780163, 2003566595, 2004090881, 2004353027, 2004877315,
            2006974467, 2007760899, 2008285187, 2008547331, 2013265923, 2013478915,
            2015363075, 2016149507, 2016673795, 2016886787, 2016935939, 2059403267,
            2061500419, 2062286851, 2062811139, 2063073283, 2067791875, 2069889027,
            2070675459, 2071199747, 2071461891, 2071986179, 2074083331, 2074869763,
            2075394051, 2075656195)

    values(file1)[(values(qc_file) %in% ok_qc)==FALSE] <- NA
    values(file2)[(values(qc_file) %in% ok_qc)==FALSE] <- NA
    values(file3)[(values(qc_file) %in% ok_qc)==FALSE] <- NA

    file14 <- 2.5 * ((file2 - file1) / (file2 + 6 * file1 - 7.5 * file3 + 1))
 
    file14[abs(file14)>1]<-NA

    yr <- substr(b01_files[i],31,32)
    y1 <- substr(y,3,4)
    DOY_layer <- paste('DOY',substr(b01_files[i],33,35),sep="")
    if (yr<y1){
      DOY_layer <- paste('DOY',as.character(-(as.numeric(substr(b01_files[i],33,35))-364)),sep = ".")
    }else if (yr==y1){
      if(substr(DOY_layer,4,5)=="00"){
        DOY_layer <- paste('DOY',substr(DOY_layer,6,6),sep="")
      }else if(substr(DOY_layer,4,4)=="0"){
        DOY_layer <- paste('DOY',substr(DOY_layer,5,6),sep="")
      }
    }else if (yr>y1){
      DOY_layer <- paste('DOY',as.character(365+as.numeric(substr(DOY_layer,6,6))),sep="")
    }
    
    EVI_raster_V061[[DOY_layer]]<-file14
    print(paste0("Load MODIS data: ",round(i/length(b01_files)*100,1),"%"))
  }
  
  rm(b01_file,b02_file,b03_file,file1,file2,file3,file14,qc_file)
  
  Inter_EVI_raster_V061 <-raster()
  #set the number of columns, rows, and extent
  Inter_EVI_raster_V061 <- brick(ncol=ncol(lc.rt), nrow=nrow(lc.rt), xmn=minlon, xmx=maxlon, ymn=minlat, ymx=maxlat,nl=380)
  res(Inter_EVI_raster_V061)
  #check the number of cells is 344448 for the GTA
  ncell(Inter_EVI_raster_V061)
  names(Inter_EVI_raster_V061)<-paste('DOY',-3:376,sep="") #name the layers

  
  EVI_values_V061<-values(EVI_raster_V061)
  EVI_length_V061<-length(values(EVI_raster_V061$DOY1))
  Inter_EVI_values_V061<-values(Inter_EVI_raster_V061)
  
  DOY<-c(-3:376)
  
  print("Interpolate MODIS Data")
  interpolate <- function(evi){
    if(sum(evi,na.rm=TRUE)>0){
      pix<-data.frame(DOY,evi)
      names(pix)<-c('DOY','EVI')
      spl<-with(pix[!is.na(pix$EVI),],smooth.spline(DOY,EVI,spar = .25))
      evi_inter<-predict(spl,DOY)$y
      return(evi_inter)
    }
  }
  
  library("foreach")
  
  system.time({foreach::foreach(i = 1:ncell(Inter_EVI_raster_V061)) %do% {
    EVI<-EVI_values_V061[i,]
    pix<-data.frame(DOY,EVI)
    names(pix)<-c('DOY','EVI')
    if(abs(sum(EVI,na.rm=TRUE))>0){
      spl<- with(pix[!is.na(pix$EVI),],smooth.spline(DOY,EVI, spar = .25)) #.25
      Inter_EVI_values_V061[i,]<-predict(spl, c(-3:374))$y #Change to -3:374 for non-leap year, -3:376 for leap year
    }
  }})
  values(Inter_EVI_raster_V061)<-Inter_EVI_values_V061
  #  
  
  if (as.numeric(y)%%4==0){
    #print(y)
    Inter_EVI_raster<-dropLayer(Inter_EVI_raster_V061,c(1,2,3,4,371:380)) # select only data in the current year 
  }else{
    Inter_EVI_raster<-dropLayer(Inter_EVI_raster_V061,c(1,2,3,4,370:380)) # select only data in the current year, 370 for non-leap year
  }
  
  #  #THIS WORKS!!! :D
  writeRaster(Inter_EVI_raster,filename=paste0(evi.path,"/",evi.pattern,".tif"),
              overwrite=TRUE)
}
