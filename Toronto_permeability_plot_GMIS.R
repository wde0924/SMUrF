memory.limit(size=5e5)
library("raster")
library("ggplot2")
library("shapefiles")
library("rgdal")

imported_raster_GMIS=raster('E:/Research/Impermeable_Surface_data/GMIS_HBASE/gmis_impervious_surface_percentage/17T_gmis_impervious_surface_percentage_geographic_30m.tif')
imported_raster_GMIS_err=raster('E:/Research/Impermeable_Surface_data/GMIS_HBASE/gmis_standard_error_of_impervious_surface_percentage/17T_gmis_standard_error_of_impervious_surface_percentage_geographic_30m.tif')

#imported_raster_MSMT_RF=raster('Global_Impervious_Surfaces_products/ImperviousMap_W80N50.tif')
#imported_raster_FROM_GLC=raster('FROM_GLC_impervious.tif')
#imported_raster_GHS=raster('GHS_Built_S_E2018/GHS_BUILT_S_E2018_GLOBE_R2023A_54009_10_V1_0_R4_C12.tif')
#imported_raster_GHS_C=raster('GHS_Built_C_E2018/GHS_BUILT_C_MSZ_E2018_GLOBE_R2023A_54009_10_V1_0_R4_C12.tif')
imported_raster_Toronto=raster('permeability_wgs84.tif') 
#last updated 2019 not sure what year's data it was created from
imported_raster_aci_2021=raster('E:/Research/Impermeable_Surface_data/ACI/aci_2021_on.tif') #updated yearly
imported_raster_SOLRIS=raster('SOLRIS_V3/SOLRIS_Version_3_0_LAMBERT.tif')
#Only covers data from 2000-2015

#imported_raster$permeability_wgs84
plot(imported_raster_GMIS)
plot(imported_raster_GMIS_err)
plot(imported_raster_Toronto)
plot(imported_raster_aci_2021)
plot(imported_raster_SOLRIS)

GMIS_crs = crs(imported_raster_GMIS) #same as MODIS!
Toronto_crs='+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs'
aci_crs = '+proj=aea +lat_0=40 +lon_0=-96 +lat_1=44.75 +lat_2=55.75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
SOLRIS_crs = '+proj=lcc +lat_0=0 +lon_0=-85 +lat_1=44.5 +lat_2=53.5 +x_0=930000 +y_0=6430000 +ellps=GRS80 +units=m +no_defs'
MODIS_crs = '+proj=longlat +datum=WGS84 +no_defs'

#bound_box_0 = readOGR(dsn="E:/Research/UrbanVPRM/dataverse_files/TPD/shapefiles", layer='TPD_30m_BB_4km')
bound_box_0 = readOGR(dsn="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/dataverse_files/shapefiles", layer='Entire_GTA_500m_BB')
#bound_box_0 = raster::extent(-80.9, -78.3, 42.4, 44.7) # regional extent

bound_box_Toronto = spTransform(bound_box_0, Toronto_crs)
bound_box_GMIS = spTransform(bound_box_0, GMIS_crs)

#bound_box = readOGR(dsn="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/dataverse_files/TPD/shapefiles",layer='TPD_30m_BB_4km') # shapefile data in /urbanVPRM_30m/shapefiles/
#bound_box_aci = spTransform(bound_box, aci_crs)
bound_box_SOLRIS = spTransform(bound_box_0, SOLRIS_crs)
#bound_box_3km = readOGR(dsn="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/dataverse_files/TPD/shapefiles",layer='TPD_30m_BB_2km') # shapefile data in /urbanVPRM_30m/shapefiles/
bound_box_aci = spTransform(bound_box_0, aci_crs)

Torontocrop = crop(imported_raster_Toronto,bound_box_Toronto)
#MSMT_RF_crop = crop(imported_raster_MSMT_RF,bound_box_0)
GMIScrop = crop(imported_raster_GMIS,bound_box_GMIS)
GMISerr_crop = crop(imported_raster_GMIS_err,bound_box_GMIS)
#GHScrop = crop(imported_raster_GHS,bound_box_GHS)
#GHS_Ccrop = crop(imported_raster_GHS_C,bound_box_GHS)
ACIcrop_2021 = crop(imported_raster_aci_2021,bound_box_aci)
SOLRIScrop = crop(imported_raster_SOLRIS,bound_box_SOLRIS)


#No difference between cropping first then aggregating and vice-versa
#Toronto_aggregated1<-aggregate(imported_raster_Toronto,60,fun=mean)
#Toronto_aggregatedcrop<-crop(Toronto_aggregated1,bound_box_Toronto)

Toronto_aggregated<-aggregate(Torontocrop,60,fun=mean) #aggregate to 30m resolution so it is easier to work with
Toronto_aggregated<-(1-Toronto_aggregated)*100 #change permeable area to % impermeable
#Toronto_imperm=Torontocrop
#Toronto_imperm[Torontocrop==0]=100
#Toronto_imperm[Torontocrop==1]=0

GMIS_imperm=GMIScrop
GMIS_imperm[GMIS_imperm==200]=0
GMIS_imperm[GMIS_imperm==255]=NA

GMISerr_imperm=GMISerr_crop
GMISerr_imperm[GMISerr_imperm==255]=NA

par(mar=c(3,3,3,0))
plot(GMIS_imperm, main='GMIS Impermeable Area')

par(mar=c(3,3,3,0))
plot(GMISerr_imperm, main='GMIS Impermeable Area Standard Error')

GMIS_frac_err=GMIS_imperm
GMIS_frac_err[GMIS_imperm==0]=NA
GMIS_frac_err=GMISerr_imperm/GMIS_frac_err
#GMIS_frac_err[GMIScrop==0]=NA
plot(GMIS_frac_err) # Average fractional error of non-zero GMIS data is 0.68, 
                    # median=0.61, but values are very variable

SOLRIS_imperm=SOLRIScrop
SOLRIS_imperm[SOLRIScrop==203]=100
SOLRIS_imperm[SOLRIScrop==201]=100
SOLRIS_imperm[SOLRIScrop==250]=SOLRIS_imperm[SOLRIScrop==250]*NA
SOLRIS_imperm[SOLRIS_imperm!=100 & SOLRIScrop!=250]=0
SOLRIS_imperm[is.na(SOLRIScrop)]=-999 #set US values to -999 to identify them

par(mar=c(3,3,3,0))
plot(SOLRIS_imperm, main='SOLRIS Impermeable Area')
plot(SOLRIScrop==250, main='SOLRIS Unclassified Area')
plot(SOLRIScrop==202, main='SOLRIS Previous Urban')

#create a raster
x <-raster()
#set the number of columns, rows, and extent
# 2km res x <- raster(ncol=210, nrow=170, xmn=1270260, xmx=1275360, ymn=611160, ymx=617460)
#x <- raster(ncol=115, nrow=140, xmn=1260510, xmx=1263960, ymn=419340, ymx=423540)
#GTA 500m
#x <- raster(ncol=1798, nrow=1901, xmn=1300980, xmx=1358010, ymn=529830, ymx=583770)
#TP39 4km
#x <- raster(ncol=230, nrow=279, xmn=1273200, xmx=1280100, ymn=428760, ymx=437130)
#TPD 4km
#x <- raster(ncol=230, nrow=278, xmn=1258800, xmx=1265700, ymn=417270, ymx=425610)
#Borden 4km
#x <- raster(ncol=226, nrow=278, xmn=1269420, xmx=1276200, ymn=610140, ymx=618480)
#GTA500m
x <- raster(ncol=8696, nrow=9822, xmn=1188780, xmx=1449660, ymn=390210, ymx=684870)
res(x)
#change the resolution
res(x) <- 30
res(x)
#check the number of cells is 35700 or 62828 for Borden 4km res
#For TPD 35776 for 3km or 16100 for 2km
#For entire GTA: 85412112
ncell(x)

# set the coordinate reference system (CRS) (define the projection)
projection(x) <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=44.75 +lat_2=55.75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
# give x the same values as ACIcrop
values(x)<-values(ACIcrop_2021)
x

ACI_imperm_2021=x
ACI_imperm_2021[x==34]=100
ACI_imperm_2021[x==35]=100 #Greenhouses, not sure if this should be impermeable or not...
ACI_imperm_2021[ACI_imperm_2021!=100]=0
rm(x)

par(mar=c(3,3,3,0))
plot(ACI_imperm_2021, main='ACI Urban and Greenhouses')
plot(ACIcrop_2021==35, main='ACI Greenhouses') 
#(barren includes rock,mines,rubble, and natural non-vegetated surfaces)')

aci_crs2="+proj=aea +datum=WGS84 +units=m +no_defs"

#SOLRIS_imperm_AEA = projectRaster(SOLRIS_imperm, crs=aci_crs, method='ngb')
#SOLRIS_imperm_crop = mask(SOLRIS_imperm_AEA, bound_box_aci, method='ngb')
#ACI_imperm_lcc = projectRaster(ACI_imperm, crs=SOLRIS_crs, method='ngb')
ACI_imperm_WGS84_2021 = projectRaster(ACI_imperm_2021, crs=MODIS_crs, method='ngb')

#rm(ACI_imperm_lcc)

par(mar=c(3,3,3,0))
##plot(SOLRIS_imperm_AEA, main='SOLRIS Impermeable Surface (proj=AEA)')
#plot(ACI_imperm_lcc, main='ACI Impermeable Surface (proj=lcc)')
plot(ACI_imperm_WGS84_2021, main='ACI 2021 Impermeable Surface (proj=WGS84)')

##ACI_imperm_AEA = projectRaster(ACI_imperm, crs=aci_crs, method='ngb')
##ACI_imperm_crop = mask(ACI_imperm_AEA, bound_box_aci, method='ngb')
##par(mar=c(3,3,3,0))
##plot(ACI_imperm_crop, main='ACI Impermeable Surface (proj=AEA)')

##ACI_resample<-resample(ACI_imperm,ACIcrop, method='ngb')
##plot(ACI_resample,main='ACI Impermeable Surface (proj=AEA)')

##SOLRIS_resample<-resample(SOLRIS_imperm_AEA,ACIcrop, method='ngb')
##plot(SOLRIS_resample, main='SOLRIS impermeable Surface (proj=AEA)')
##SOLRIS_imperm_crop = crop(SOLRIS_resample, bound_box_aci, method='ngb')
##plot(SOLRIS_imperm_crop,main='SOLRIS Impermeable Surface (proj=AEA)')

#SOLRIS_crop_WGS84 = projectRaster(SOLRIScrop, crs=MODIS_crs, method='ngb')
#SOLRIS_crop_WGS84 <- resample(SOLRIS_crop_WGS84,GMIScrop,method='ngb')
#plot(SOLRIS_crop_WGS84, main='SOLRIS (proj=WGS84)')

SOLRIS_imperm_WGS84 = projectRaster(SOLRIS_imperm, crs=MODIS_crs, method='ngb')
#If using GMIS uncomment this line:
SOLRIS_resample_WGS84_GMIS <- resample(SOLRIS_imperm_WGS84,GMIScrop,method='ngb')
#If not using GMIS uncomment this line:
#SOLRIS_resample_WGS84 <- resample(SOLRIS_imperm_WGS84,ACI_imperm_WGS84,method='ngb')
plot(SOLRIS_resample_WGS84_GMIS, main='SOLRIS impermeable Surface (proj=WGS84)')

#ACI_resample<-resample(ACI_imperm_lcc,SOLRIScrop, method='ngb')
#ACI_resample <- resample(ACI_imperm_lcc,SOLRIS_imperm,method='ngb')
#plot(ACI_resample, main='ACI impermeable Surface (proj=lcc)')

#ACI does not have any NA values, only 0's

ACI_resample_WGS84_GMIS_2021 <- resample(ACI_imperm_WGS84_2021,GMIScrop,method='ngb')
plot(ACI_resample_WGS84_GMIS_2021, main='ACI 2021 impermeable Surface (proj=WGS84)')

#SOLRIS_imperm[SOLRIScrop==250]=ACI_resample[SOLRIScrop==250]
#plot(SOLRIS_imperm, main='SOLRIS Impermeable Area')

#Toronto_imperm_lcc = projectRaster(Toronto_aggregated, crs=SOLRIS_crs, method='ngb')
#Toronto_resample<-resample(Toronto_imperm_lcc,SOLRIScrop,method='ngb')
#plot(Toronto_resample,main='Toronto Impermeable Area (proj=lcc)')

Toronto_imperm_wgs84 = projectRaster(Toronto_aggregated, crs=MODIS_crs, method='ngb')
#If using GMIS:
Toronto_resample_wgs84_GMIS<-resample(Toronto_imperm_wgs84,GMIScrop,method='ngb')
#If not using GMIS:
#Toronto_resample_wgs84<-resample(Toronto_imperm_wgs84,ACI_imperm_WGS84,method='ngb')
plot(Toronto_resample_wgs84_GMIS,main='Toronto Impermeable Area (proj=WGS84)')

#Anywhere considered urban by ACI is replaced by SOLRIS urban pervious or urban impervious
ACI_SOLRIS_WGS84_GMIS_2021 <- ACI_resample_WGS84_GMIS_2021 #ACI_imperm_WGS84
ACI_SOLRIS_WGS84_GMIS_2021[SOLRIS_resample_WGS84_GMIS==-999]=NA
ACI_SOLRIS_WGS84_GMIS_2021[ACI_resample_WGS84_GMIS_2021==100 & is.na(SOLRIS_resample_WGS84_GMIS)==FALSE]=SOLRIS_resample_WGS84_GMIS[ACI_resample_WGS84_GMIS_2021==100 & is.na(SOLRIS_resample_WGS84_GMIS)==FALSE]
ACI_SOLRIS_WGS84_GMIS_2021[ACI_SOLRIS_WGS84_GMIS_2021==-999]=NA
plot(ACI_SOLRIS_WGS84_GMIS_2021,main='2021 ACI & SOLRIS ISA')

#ACI_resample[ACI_resample==100]<-75
#plot(ACI_resample,main='ACI & SOLRIS Impermeable Area with 75% impervious')

#Toronto_imperm_WGS84 <- projectRaster(Toronto_aggregated, crs=MODIS_crs)
#GMIS_resample_WGS84 <- resample(GMIS_imperm, SOLRIS_resample_WGS84_GMIS, method='ngb')
#plot(Toronto_resample_wgs84_GMIS, main='Toronto Impermeable Percentage (proj=WGS84)')
#plot(GMIS_resample_WGS84, main='GMIS Percentage (proj=WGS84)')
#plot((Toronto_resample_wgs84-GMIS_resample_WGS84),main='Toronto impermeable - GMIS %')
##On average Toronto impermeable surface pixels are 6.59% more impervious
## Median difference is 3.59%

med_Toronto_SOLRIS_2021<-median(Toronto_resample_wgs84_GMIS[ACI_SOLRIS_WGS84_GMIS_2021==100],na.rm=TRUE) #63.97% #sd=29.27%
mean_Toronto_SOLRIS_2021<-mean(Toronto_resample_wgs84_GMIS[ACI_SOLRIS_WGS84_GMIS_2021==100],na.rm=TRUE) #60.52%

med_GMIS_SOLRIS_2021<-median(GMIS_imperm[ACI_SOLRIS_WGS84_GMIS_2021==100],na.rm=TRUE) #31% #sd=32.09%
mean_GMIS_SOLRIS_2021<-mean(GMIS_imperm[ACI_SOLRIS_WGS84_GMIS_2021==100],na.rm=TRUE) #35.13% 

mean_GMIS_ACI<-mean(GMIS_imperm[ACI_resample_WGS84_GMIS_2021==100],na.rm=TRUE) #30.18% #sd=31.82%
med_GMIS_ACI<-median(GMIS_imperm[ACI_resample_WGS84_GMIS_2021==100],na.rm=TRUE) #22%

Toronto_ACI_SOLRIS_2021<-ACI_SOLRIS_WGS84_2021
Toronto_ACI_SOLRIS_2021[ACI_SOLRIS_WGS84_2021==100 & is.na(Toronto_resample_wgs84)]=med_Toronto_SOLRIS_2021
Toronto_ACI_SOLRIS_2021[is.na(Toronto_resample_wgs84)==FALSE]<-Toronto_resample_wgs84[is.na(Toronto_resample_wgs84)==FALSE]
plot(Toronto_ACI_SOLRIS_2021,main='2021 Toronto, ACI, & SOLRIS Impervious (proj=WGS84)')

GMIS_ACI_SOLRIS_2021<-ACI_SOLRIS_WGS84_GMIS_2021
GMIS_ACI_SOLRIS_2021[ACI_SOLRIS_WGS84_GMIS_2021==100 & GMIS_imperm==0]=med_GMIS_SOLRIS_2021
GMIS_ACI_SOLRIS_2021[ACI_SOLRIS_WGS84_GMIS_2021==100 & GMIS_imperm>0]=GMIS_imperm[ACI_SOLRIS_WGS84_GMIS_2021==100 & GMIS_imperm>0]
GMIS_ACI_SOLRIS_2021[is.na(ACI_SOLRIS_WGS84_GMIS_2021)]<-GMIS_imperm[is.na(ACI_SOLRIS_WGS84_GMIS_2021)]
GMIS_ACI_SOLRIS_2021[GMIS_ACI_SOLRIS_2021==-999]=NA
plot(GMIS_ACI_SOLRIS_2021,main='GMIS, 2021 ACI, & SOLRIS Impervious (proj=WGS84)')

GMIS_ACI<-ACI_resample_WGS84
GMIS_ACI[ACI_resample_WGS84==100 & GMIS_imperm==0]=med_GMIS_ACI
GMIS_ACI[ACI_resample_WGS84==100 & GMIS_imperm>0]=GMIS_imperm[ACI_resample_WGS84==100 & GMIS_imperm>0]
GMIS_ACI[is.na(ACI_SOLRIS_WGS84)]<-GMIS_imperm[is.na(ACI_SOLRIS_WGS84)]
plot(GMIS_ACI,main='GMIS & ACI Impervious (proj=WGS84)')


GMIS_Toronto_ACI_SOLRIS_2021<-GMIS_ACI_SOLRIS_2021
GMIS_Toronto_ACI_SOLRIS_2021[is.na(Toronto_resample_wgs84_GMIS)==FALSE]<-Toronto_resample_wgs84_GMIS[is.na(Toronto_resample_wgs84_GMIS)==FALSE]
plot(GMIS_Toronto_ACI_SOLRIS_2021,main='GMIS-Toronto-ACI2021-SOLRIS ISA')
#Average difference between GMIS_ACI & GMIS_ACI_SOLRIS is 0.18% with sd:3.60%

#Toronto_imperm_moll_mean = projectRaster(Toronto_aggregated, crs=GHS_crs)
#GHS_resample<-resample(GHScrop,Toronto_imperm_moll,method='ngb')
#GHS_C_resample<-resample(GHS_Ccrop,Toronto_imperm_moll,method='ngb')
#plot(Toronto_imperm_moll,main='Toronto Impermeable Area (proj=moll)')
#plot(GHS_resample,main='GHS Built-up fraction (proj=moll)')

#GHS_imperm_merc = projectRaster(GHScrop, crs=Toronto_crs, method='ngb')
#GHS_resample_merc<-resample(GHS_imperm_merc,Toronto_aggregated,method='ngb')
#plot(Toronto_aggregated,main='Toronto Impermeable Area (proj=merc)')
#plot(GHS_resample_merc,main='GHS Built-up fraction (proj=merc)')

#LC = raster('E:/Research/UrbanVPRM/dataverse_files/Borden_500m/LandCover/MODIS_LC_Borden_500m.tif') # Land cover data in /urbanVPRM_30m/driver_data/lc_isa/
LC = raster('C:/Users/kitty/Documents/Research/SIF/SMUrF/data/MCD12Q1/MCD12Q1.061_LC_Type1_doy2021001_aid0001.tif')
crs(LC)<-MODIS_crs

#Optional (only needed if using Toronto dataset without GMIS)

##LC_llc<-projectRaster(LC,crs=SOLRIS_crs,method='ngb')
#LC_resample<-resample(LC,SOLRIS_imperm_WGS84,method='ngb')

#Toronto_WGS84_crop<-crop(Toronto_resample_wgs84,LC_resample)
#ACI_WGS84_crop<-crop(ACI_resample,LC_resample)
#ACI_WGS84_crop[values(is.na(Toronto_WGS84_crop))] <- NA
#plot(ACI_WGS84_crop)
#plot(Toronto_WGS84_crop-ACI_WGS84_crop)

Toronto_vals<-values(Toronto_resample_wgs84)#values(Toronto_WGS84_crop)
ACI_vals<-values(ACI_resample_WGS84) #values(ACI_lcc_crop)
plot(Toronto_vals[ACI_vals==100],ACI_vals[ACI_vals==100])
#check the median value of ISA from the Toronto dataset when the ISA from ACI is 100
median(Toronto_vals[ACI_vals==100],na.rm=TRUE) #63.47 (mean=59.85, sd=29.55)
median(Toronto_vals[ACI_vals==0],na.rm=TRUE) #1.53 (mean=29.55, sd=30.9)

#end of optional

#Set ACI 100% ISA values to the median Toronto values
ACI_resample[ACI_resample==100]<-63
plot(ACI_resample,main='ACI & SOLRIS Impermeable Area with 63% impervious')

plot(Toronto_resample-ACI_resample,main='Toronto ISA - ACI & SOLRIS ISA with 63% impervious')
#The median difference in values is 1.17 (mean=1.25, sd=31.04)
#If use 75% instead get median difference of -8.81 (mean=-14.29,SD=34.56)
#If use 60% instead get median difference of 0.61 (mean=-1.21, SD=33.04)



#Toronto_resample_mod<-projectRaster(Toronto_aggregated, crs=crs(LC))
##Toronto_resample_mod<-crop(Toronto_resample_mod,LC)
LC_crop <- crop(LC,bound_box_0)
#Toronto_resample_mod <- crop(Toronto_resample_mod,LC_crop)
#Toronto_resample_mod <- resample(Toronto_resample_mod,LC_crop)
GMIS_resample_mod <- resample(GMIS_imperm,LC_crop)
GMIS_resample_mod[GMIS_resample_mod<0]<-NA
GMISerr_resample_mod <- resample(GMISerr_imperm,LC_crop)
GMISerr_resample_mod[GMISerr_resample_mod<0] <- NA
plot(GMIS_resample_mod, main='Aggregated GMIS ISA %')
plot(GMISerr_resample_mod, main='Aggregated GMIS ISA % Error')

Toronto_resample_mod <- resample(Toronto_resample_wgs84,LC_crop)
plot(Toronto_resample_mod, main='Aggregated City of Toronto Impermeable Surface area %')

Toronto_ACI_SOLRIS_mod_2021 <- resample(Toronto_ACI_SOLRIS_2021,LC_crop)
plot(Toronto_ACI_SOLRIS_mod_2021, main='Aggregated 2021 Toronto-ACI-SOLRIS Impermeable Surface area %')


GMIS_ACI_SOLRIS_resample_mod_2021 <- resample(GMIS_ACI_SOLRIS_2021,LC_crop)
GMIS_ACI_SOLRIS_resample_mod_2021[GMIS_ACI_SOLRIS_resample_mod_2021<0]<-NA
plot(GMIS_ACI_SOLRIS_resample_mod_2021, main='Aggregated GMIS-2021ACI-SOLRIS ISA %')

GMIS_ACI_resample_mod <- resample(GMIS_ACI, LC_crop)
GMIS_ACI_resample_mod[GMIS_ACI_resample_mod<0]<-NA
plot(GMIS_ACI_resample_mod, main='Aggregated GMIS-ACI Impermeable Surface area %')

GMIS_Toronto_ACI_SOLRIS_resample_mod_2021 <- resample(GMIS_Toronto_ACI_SOLRIS_2021, LC_crop)
GMIS_Toronto_ACI_SOLRIS_resample_mod_2021[GMIS_Toronto_ACI_SOLRIS_resample_mod_2021<0]<-0
plot(GMIS_Toronto_ACI_SOLRIS_resample_mod_2021, main='Aggregated GMIS-Toronto-2021ACI-SOLRIS ISA %')

#SOLRIS_imperm[SOLRIS_imperm==100]<-63
#SOLRIS_resample_mod<-projectRaster(SOLRIS_imperm, crs=crs(LC))
#SOLRIS_resample_mod<-crop(SOLRIS_resample_mod,LC_crop)
#SOLRIS_resample_mod <- resample(SOLRIS_resample_mod,LC_crop)
#plot(SOLRIS_resample_mod, main='Aggregated SOLRIS Impermeable Surface area %')

#ACI_imperm[ACI_imperm==100]<-63
#LC_reproj <- projectRaster(LC_crop,crs=crs(ACI_imperm))
#ACI_imperm_crop  <- resample(ACI_imperm, LC_reproj)
#ACI_resample_mod<-projectRaster(ACI_imperm_crop, crs=crs(LC))
#ACI_resample_mod<-crop(ACI_resample_mod,LC_crop)
#ACI_resample_mod <- resample(ACI_resample_mod,LC_crop)
#plot(ACI_resample_mod, main='Aggregated ACI Impermeable Surface area %')


#ACI_SOLRIS_resample_mod<-projectRaster(ACI_resample, crs=crs(LC))
#ACI_SOLRIS_resample_mod<-crop(ACI_SOLRIS_resample_mod,LC_crop)
#ACI_SOLRIS_resample_mod <- resample(ACI_SOLRIS_resample_mod,LC_crop)
#plot(ACI_SOLRIS_resample_mod, main='Aggregated ACI & SOLRIS(63%) Impermeable Surface area %')
##median difference between aggregated Toronto and ACI & SOLRIS ISA is 0.49 (mean=3.21, sd=14.09)
##using 75% median difference is -14.53 (mean=-11.68, sd=16.87)

##When Toronto impermeable surface is not available replace it with ACI values (impermeable )
all_resample_mod<-GMIS_ACI_SOLRIS_resample_mod
all_resample_mod[!is.na(Toronto_resample_mod)]=Toronto_resample_mod[!is.na(Toronto_resample_mod)]
all_resample_mod[all_resample_mod<0]<-NA
#all_resample<-Toronto_resample_mod
#all_resample[is.na(Toronto_resample_mod) & is.na(SOLRIS_resample_mod)==FALSE]<-SOLRIS_resample_mod[is.na(Toronto_resample_mod) & is.na(SOLRIS_resample_mod)==FALSE]
#all_resample[is.na(Toronto_resample_mod) & is.na(SOLRIS_resample_mod)] <- ACI_resample_mod[is.na(Toronto_resample_mod) & is.na(SOLRIS_resample_mod)]
#all_resample[is.na(Toronto_resample_mod)]<-ACI_SOLRIS_resample_mod[is.na(Toronto_resample_mod)]
plot(all_resample_mod, main='Combined Impermeable Surface Area %')
#I am not happy with the area surrounding Toronto, 
# it appears too impervious compared to the city

writeRaster(GMIS_resample_mod,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/GMIS_impervious_GTA.tif",
            overwrite=TRUE)
writeRaster(GMISerr_resample_mod,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/GMISerr_impervious_GTA.tif",
            overwrite=TRUE)
writeRaster(Toronto_ACI_SOLRIS_mod_2018,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/all_30m_aggregated_2018_impervious_63_GTA.tif",
            overwrite=TRUE)
writeRaster(Toronto_ACI_SOLRIS_mod,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/all_30m_aggregated_impervious_63_GTA.tif",
            overwrite=TRUE)
writeRaster(Toronto_ACI_SOLRIS_mod_2020,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/all_30m_aggregated_2020_impervious_63_GTA.tif",
            overwrite=TRUE)
writeRaster(GMIS_ACI_SOLRIS_resample_mod_2018,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/GMIS_ACI_SOLRIS_2018_impervious_GTA.tif",
            overwrite=TRUE)
writeRaster(GMIS_Toronto_ACI_SOLRIS_resample_mod_2021,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/GMIS_Toronto_ACI_SOLRIS_2021_impervious_GTA.tif",
            overwrite=TRUE)

writeRaster(all_resample_mod,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/Toronto_GMIS_ACI_SOLRIS_impervious_GTA.tif",
            overwrite=TRUE)
writeRaster(GMIS_ACI_resample_mod,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/GMIS_ACI_impervious_GTA.tif",
            overwrite=TRUE)
#writeRaster(SOLRIS_resample_mod,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/SOLRIS_aggregated_impervious_GTA.tif",
#            overwrite=TRUE)

writeRaster(GMIS_ACI_SOLRIS,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/GMIS_ACI_SOLRIS_30m_impervious_GTA.tif",
            overwrite=TRUE)

GMIS_Toronto_ACI_SOLRIS <- raster("C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/GMIS_Toronto_ACI_SOLRIS_2018_impervious_GTA.tif")
plot(GMIS_Toronto_ACI_SOLRIS,main='Toronto, GMIS, SOLRIS, & ACI Impermeable Surface %')


#Look at ISA without GMIS data: 
all_resample<-raster("C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/all_30m_aggregated_2018_impervious_63_GTA.tif")
plot(all_resample,main='Toronto, SOLRIS, & ACI Impermeable Surface %')

#all_aggregated<-aggregate(all_resample,240/20,fun=mean) #aggregate to 0.05 degree resolution so it is easier to work with
#plot(all_aggregated,main='All aggregated ISA, 0.05^o resolution')

#writeRaster(all_aggregated,filename="C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface/all_aggregated_impervious_63_entire_GTA_CSIF_res.tif",
#            overwrite=TRUE)
