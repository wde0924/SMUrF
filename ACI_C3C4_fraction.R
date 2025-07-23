library("raster")
library("ggplot2")
library("shapefiles")
library("rgdal")

homedir <- 'C:/Users/kitty/Documents/Research/SIF'
smurf_wd <- file.path(homedir, 'SMUrF'); setwd(smurf_wd)
source('r/dependencies.r')

imported_raster_aci=raster('E:/Research/Impermeable_Surface_data/ACI/aci_2018_qc_v3.tif')
imported_raster_aci2=raster('E:/Research/Impermeable_Surface_data/ACI/aci_2018_on_v3.tif')

plot(imported_raster_aci)
plot(imported_raster_aci2)

aci_crs = '+proj=aea +lat_0=40 +lon_0=-96 +lat_1=44.75 +lat_2=55.75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

indx<-2
reg.name <- c('westernCONUS', 'easternCONUS',     'westernEurope', 
              'easternChina', 'easternAustralia', 'easternAsia', 
              'southAmerica', 'centralAfrica')[indx]

# Southern Ontario: -80.9, -78.3, 42.4, 44.7

minlon <- c(-125, -76.5,  -11, 100,  130, 125, -65, -10)[indx]
maxlon <- c( -95, -72.7,   20, 125,  155, 150, -40,  20)[indx]
minlat <- c(  25,  44.5,   35,  20,  -40,  30, -40, -10)[indx]
maxlat <- c(  50,  46.4,   60,  50,  -10,  55, -10,  15)[indx]
reg.ext<-extent(minlon, maxlon, minlat, maxlat)

lc.path    <- file.path(smurf_wd, 'data/MCD12Q1/Montreal_Ottawa')
lc.pattern <- 'MCD12Q1.061_LC_Type1'

# indicate the latest year available of MCD12Q1
# if no data beyond 2022, use 2022 LC for 2023 and beyond
lc.max.yr <- 2018  
lc.res    <- 1/24

lc.rt <- prep.mcd12(lc.path, lc.pattern, 2018, lc.max.yr, reg.name, reg.ext)

p <- as(reg.ext, 'SpatialPolygons')
crs(p) <- "+proj=longlat +datum=WGS84 +no_defs"
bound_box_aci<-spTransform(p,aci_crs)

ACIcrop <-crop(imported_raster_aci,bound_box_aci)
ACIcrop2 <-crop(imported_raster_aci2,bound_box_aci)


merged_raster <- mosaic(ACIcrop,ACIcrop2,fun=max)

#create a raster
x <-raster()
#set the number of columns, rows, and extent
# 2km res x <- raster(ncol=210, nrow=170, xmn=1270260, xmx=1275360, ymn=611160, ymx=617460)
#x <- raster(ncol=8696, nrow=9822, xmn=1188780, xmx=1449660, ymn=390210, ymx=684870) #For S. Ontario
x <- raster(ncol=ncol(merged_raster), nrow=nrow(merged_raster), xmn=xmin(merged_raster), xmx=xmax(merged_raster), ymn=ymin(merged_raster), ymx=ymax(merged_raster)) #For Montreal/Ottawa
res(x)
## check the number of cells is 85412112 for S. Ontario
## check the number of cells is 100342422 for Montreal/Ottawa region
ncell(x) == ncell(merged_raster)

# set the coordinate reference system (CRS) (define the projection)
projection(x) <- aci_crs
## give x the same values as ACIcrop/merged raster
values(x)<-values(merged_raster)
#change values outside the domain to -100 and NA over water
x[x==0]<- -100  
x[x==20]<- 0


ACI_LC=x
rm(x, ACIcrop, ACIcrop2, merged_raster)
#ACI_LC[x==20]=0 #Water -> 0
#ACI_LC[x==30]=0 #Exposed/Baren -> 0

#ACI_LC[x==34]=0 #?Urban/Developed -> 0
#ACI_LC[x==35]=0 #?Greenhouses -> 0?

#ACI_LC[x==50]=0 #Shrubland -> 0 # Shrubs are C3

#ACI_LC[x==80]=0 #?Wetland -> 0?
                # It seems most above ground (not underwater) wetland species
                # including mosses are C3 
# Borden are classified as swamps MAY NEED TO REVISE
# Not sure if this works for TP since there are both swamps and marshes...
# https://geohub.lio.gov.on.ca/datasets/mnrf::wetlands/explore?filters=eyJXRVRMQU5EX1RZUEUiOlsiU3dhbXAiLCJGZW4iLCJCb2ciLCJVbmtub3duIiwiTWFyc2giLCJPcGVuIFdhdGVyIl19&location=44.296891%2C-79.932475%2C13.35

#ACI_LC[x==85]=0 #?Peatland (wetland comercially harvested for peat) -> 0?
                # Typically form due to Sphagnum mosses (bogs and fens rare in S.ON)
                # Can contain herbs, ferns, grasses, sedges, shrubs and trees
#ACI_LC[ACI_LC<100]=0

#ACI_LC[ACI_LC==110]=0.95 #?Grassland -> 1? (95% of S.Ont grasses seem to be C4 remainder 
                 # is mostly Canada Rye (C3) and a few mixed species)
#only used for crops in SMUrF

#ACI_LC[ACI_LC==120]=0 #?Agriculture (undifferentiated) -> 0? (look into average % of C4 plants in crops?)

#ACI_LC[ACI_LC==122]=0 # Passture/Forages -> 0? Look into what kind of grasses are typically used
                 # Pastures are usually a mixture of grasses and legumes (clover)
                 # Most common grasses according to OMAFRA are timothy (C3), 
                 # brome (C3), orchard (C3), reed canary (C3), fescues (C3),
                 # meadow foxtail (C3), bluegrass (C3), ryegrass (C3), 
                 # Redtop (C3), and Quack (C3). Legumes are C3

#ACI_LC[ACI_LC==130]=0 #Too wet to seed -> 0
#ACI_LC[ACI_LC==131]=0 #Fallow -> 0

#ACI_LC[ACI_LC==132]=0 #?Cereals -> 0? (Most are C3 (wheat, oats, rye, rice, barley),
                 #  but corn and types of millet are C4 check average fraction)
#ACI_LC[ACI_LC==133]=0 #Barley is C3 -> 0
#ACI_LC[ACI_LC==134]=0 #?Other grains -> 0?
ACI_LC[ACI_LC==135]=1 #Millet is C4 -> 1
#ACI_LC[ACI_LC==136]=0 #Oat is C3 -> 0
#ACI_LC[ACI_LC==137]=0 #Rye is C3 -> 0
#ACI_LC[ACI_LC==138]=0 #Spelt is C3 -> 0
#ACI_LC[ACI_LC==139]=0 #Triticale is C3 -> 0 (I think)?
#ACI_LC[ACI_LC==140]=0 #Wheat is C3 -> 0
ACI_LC[ACI_LC==141]=1 #Switchgrass is C4 -> 1
ACI_LC[ACI_LC==142]=1 #Sorghum is C4 -> 1
#ACI_LC[ACI_LC==143]=0 # Quinoa is C3? -> 0 THERE SEEMS TO BE CONFLICTING INFO ON THIS
                 # MAYBE SOME VARIATIONS ARE C3 AND SOME C4? 
                 # I think the quinoa grown for human consumption is C3? 
                 # It seems to be debated though...
#ACI_LC[ACI_LC==145]=0 #Winter Wheat is C3 -> 0
#ACI_LC[ACI_LC==146]=0 #Spring Wheat is C3 -> 0
ACI_LC[ACI_LC==147]=1 #Corn is C4 -> 1
#ACI_LC[ACI_LC==148]=0 #Tobacco is C3 -> 0
#ACI_LC[ACI_LC==149]=0 #Ginseng is C3 -> 0

#ACI_LC[ACI_LC==150]=0 #main oilseeds (sunflower,canola, soy, flax) are C3 ->0
#ACI_LC[ACI_LC==151]=0 #Borage is C3 -> 0 (I think)
#ACI_LC[ACI_LC==152]=0 #Camelina is C3 -> 0
#ACI_LC[ACI_LC==153]=0 #Canola is C3 -> 0
#ACI_LC[ACI_LC==154]=0 #Flaxseed is C3 -> 0
#ACI_LC[ACI_LC==155]=0 #Mustard is C3 -> 0
#ACI_LC[ACI_LC==156]=0 #Safflower is C3 -> 0
#ACI_LC[ACI_LC==157]=0 #Sunflower is C3 -> 0
#ACI_LC[ACI_LC==158]=0 #Soybean is C3 -> 0

#ACI_LC[ACI_LC==160]=0 #Most Pulses are C3 -> 0
#ACI_LC[ACI_LC==161]=0 #Other Pulses is C3 -> 0
#ACI_LC[ACI_LC==162]=0 #Peas are C3 -> 0
#ACI_LC[ACI_LC==163]=0 #Chickpeas are C3 -> 0
#ACI_LC[ACI_LC==167]=0 #Beans are C3 -> 0
#ACI_LC[ACI_LC==168]=0 #Fababeans (fava beans) are C3 -> 0
#ACI_LC[ACI_LC==174]=0 #Lentils are C3 -> 0

#ACI_LC[ACI_LC==175]=0 #Most Vegetables are C3 -> 0
#ACI_LC[ACI_LC==176]=0 #Tomatoes are C3 -> 0
#ACI_LC[ACI_LC==177]=0 #Potatoes are C3 -> 0
#ACI_LC[ACI_LC==178]=0 #Sugarbeet is C3 -> 0
#ACI_LC[ACI_LC==179]=0 #Other Vegetables -> 0

#ACI_LC[ACI_LC==180]=0 #Most Fruit is C3 -> 0
#ACI_LC[ACI_LC==181]=0 #All Berries are C3 -> 0
#ACI_LC[ACI_LC==182]=0 #Blueberry is C3 -> 0
#ACI_LC[ACI_LC==183]=0 #Cranberry is C3 -> 0
#ACI_LC[ACI_LC==185]=0 #Other Berries are C3 -> 0
#ACI_LC[ACI_LC==188]=0 #Orchards (all trees are C3) are C3 -> 0
#ACI_LC[ACI_LC==189]=0 #Other Fruit -> 0
#ACI_LC[ACI_LC==190]=0 #Vinyards (grapes) are C3 -> 0
#ACI_LC[ACI_LC==191]=0 #Hops are C3 -> 0

#ACI_LC[ACI_LC==192]=0 #SOD is C3 (assuming cool-season grasses in Toronto)-> 0
#ACI_LC[ACI_LC==193]=0 # I AM NOT SURE ABOUT HERBS, 
                 # I THINK MOST EDDIBLE ONES ARE C3 BUT I AM NOT SURE
#ACI_LC[ACI_LC==194]=0 #Most plants are C3 (similar for nurseries) ->0
#ACI_LC[ACI_LC==195]=0 # I DON'T THINK BUCKWHEAT IS C4 (it belongs to the polygonaceae 
                 # family but the only genus in this family that is C4 is 
                 # Calligonum but Buckwheat is of the genus Fagopyrum)
#ACI_LC[ACI_LC==196]=0 #Canaryseed is a BOP grass (all of which are C3) ->0
#ACI_LC[ACI_LC==197]=0 # Hemp is C3 -> 0
#ACI_LC[ACI_LC==199]=0 # Most Crops (asside from Corn and Millet which are already 
                 # accounted for) are C3, Other Crops -> 0

#ACI_LC[ACI_LC>199]=0 # All tree species are C3 thus Forests -> 0

ACI_LC[ACI_LC>2]=0

#par(mar=c(3,3,3,0))
#plot(ACI_LC, main='ACI C4 Plants')

ACI_C4_proj<-projectRaster(ACI_LC, crs='+proj=longlat +datum=WGS84 +no_defs', method='ngb')
#par(mar=c(3,3,3,0))
#plot(ACI_C4_proj, main='ACI C4 fraction (proj=WGS84)')

ACI_C4_crop<-crop(ACI_C4_proj, p, method='ngb')
#par(mar=c(3,3,3,0))
#plot(ACI_C4_crop, main='ACI C4 fraction cropped')


ACI_resample<-resample(ACI_C4_crop,lc.rt)
#plot(ACI_resample,main='ACI C4 fraction')


#If there is no data
ratio.fn <- file.path(smurf_wd, 'data/C4_relative_fraction.tif')
SPAM_C4 <- raster(ratio.fn)
SPAM_C4_crop <- crop(SPAM_C4,p,method='ngb')
SPAM_C4_resample <- resample(SPAM_C4_crop,lc.rt,method='ngb')

##smurf_C4<-raster("E:/Research/SMUrF/output2018_500m_CSIF_to_TROPOMI_CSIF_ALL_converted_slps_temp_impervious_R_8day/easternCONUS/C4_ratio_easternCONUS.tif")



SPAM_US_crops <- SPAM_C4_resample
SPAM_US_crops[(ACI_resample>=0)]<-NA #Isolate only the areas in the US (where ACI does not exist)
SPAM_US_crops[(lc.rt!=12) & (lc.rt!=13)]<-NA #isolate only areas defined by MODIS as crops

# Calculate the average C3:C4 ratio in the region (in the US only since this is
# where the filling occurs)
SPAM_US_mean_C4_ratio <- mean(values(SPAM_US_crops),na.rm=T)


ACI_resample[ACI_resample< 0]=SPAM_C4_resample[ACI_resample< 0]
ACI_resample[is.na(ACI_resample)]=SPAM_C4_resample[is.na(ACI_resample)]


## If there are areas over croplands that are NA (due to coarser resolution of SPAM)
## set them to the average C3:C4 ratio in the region (in the US)
ACI_resample_filled <- ACI_resample
ACI_resample_filled[is.na(ACI_resample) & ((lc.rt==12) | (lc.rt==14))] <- SPAM_US_mean_C4_ratio
plot(ACI_resample_filled,main='ACI & SPAM C4 fraction')


writeRaster(ACI_resample_filled,filename="C:/Users/kitty/Documents/Research/SIF/SMUrF/data/ACI_C4_fraction_Montreal_Ottawa_500m_2018.tif",
            overwrite=TRUE)
