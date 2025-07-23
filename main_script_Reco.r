#' Main script to generate ecosystem respiration
#' @author: Dien Wu, 05/20/2019
#' ---------------------------------------------------------------------------

#' @GeneralIdeas:
#' 1. Need gridded Tair, Tsoil and SIF-based GPP. 
#'    Modeled GPP can be generated from 'main_script_GPP_temporal.r'
#' 2. Need pretrained NN model derived from FLUXNET, which can be generated from 
#'    'prep_NN_train_reco.r' and 'NN_train_reco.r'

#' @InputData:
#' 0. GPP estimates based SIF, see 'main_script_GPP_temporal.r'
#' 1. Gridded 1km x1km Tmin and Tmax from Daymet 
#' 2. Gridded Tsoil from NLDAS
#' 3. Gridded sub-categories for urban settlements e.g., from NLCD
#' ---------------------------------------------------------------------------

#' @updates, by DW:
#' 06/18/2019 incorporate slurm parallel scripts to this script
#' 08/02/2019 add NN models trained by either NLDAS+daymet (US) or ERA5 (global)
#'            also NN models trained by either FLUXNET or modeled temp + GPP 
#' ---------------------------------------------------------------------------

#' @updates, by SM:
#' 09/28/2021 updated code to work at 500 m x 500 m resolution 
#'            (no longer aggregates)
#' 11/29/2021 fixed unrepresentative savannas in northern latitudes by setting 
#'            to agb-based land cover (similar to urban)
#' 11/10/2022 Adjusted R calculation to limit heterotrophic soil respiration 
#'            over impervious surfaces (following Hardimann et al. 2017)
#' 03/13/2024 Adjusted to use V6.1 of MODIS
#' 03/13/2024 Adjusted autotrophic respiration in urban areas (following 
#'            Hardimann et al. 2017)

# when using runthem.py, turn this on
#args <- commandArgs(trailingOnly = TRUE)

memory.limit(size=5e5)

#### source all functions and load all libraries
homedir <- 'C:/Users/kitty/Documents/Research/SIF'
smurf_wd <- file.path(homedir, 'SMUrF'); setwd(smurf_wd)
source('r/dependencies.r')              # source all functions


# ---------------------------------------------------------------------------
# Paths one needs to modify 
# ---------------------------------------------------------------------------
input.path  <- file.path(homedir, 'SMUrF/data')

#Output path should be the same as output for GPP
output.path <- file.path(homedir, 'SMUrF/output2018_GTA_500m_TROPOMI_CSIF_impervious_R_shore_corr_V061_8day')


## path for the updated 500m IGBP generated from main_script_GPP.r
lc.path <- file.path(smurf_wd, 'data/MCD12Q1/')
lc.pattern <- 'MCD12Q1.061_LC_Type1'

lc.max.yr<-2023

tmpdir <- 'C:/Users/kitty/AppData/Local/Temp/R' #NA

ISA.path <- 'C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/Impermeable_Surface'
ISA.pattern <- 'GMIS_Toronto_ACI_SOLRIS_2018_impervious_GTA'

isa.max.yr<-2022

# isa.tf Can be set to false for run without ISA adjustment. BUT should ALWAYS
# be set to True unless trying to recreate fig. 4 of Madsen-Colford et al. 2025
isa.tf <- TRUE
            
# ---------------------------------------------------------------------------
# Variables one needs to modify 
# ---------------------------------------------------------------------------
# name your region, needs to be consistent with that in main_script_GPPv2.r
#indx <- as.numeric(args[1])     # get region indx from python code, e.g., 1
indx <- 2
reg.name <- c('westernCONUS', 'easternCONUS',     'westernEurope', 
              'easternChina', 'easternAustralia', 'easternAsia', 
              'southAmerica', 'centralAfrica')[indx]   
              
# output path for the target region
reg.path <- file.path(output.path, reg.name)
cat(paste('Working on', reg.name, '...\n'))

# please make sure this domain is <= than the domain of MODIS land cover,
# 'minlon maxlon, minlat, laxlat' should matche the order of 'reg.name' above
# *** too large a spatial extent may lead to memory issue, DONT DO ENTIRE GLOBE
# For Southern Ontario: -80.9, -78.3, 42.4, 44.7
# For Montreal & Ottawa region: -76.2, -72.7, 44.5, 46.4
# Just greater Montreal area: -74.4, -72.9, 45.1, 46.1
# Just Canadian Capitol Region (Ottawa/Gatineau): -76.5, -75.1, 44.8, 46
#minlon <- c(-125, -74.4,  -11, 100,  130, 125, -65, -10)[indx]
#maxlon <- c( -95, -72.9,   20, 125,  155, 150, -40,  20)[indx]
#minlat <- c(  25,  45.1,   35,  20,  -40,  30, -40, -10)[indx]
#maxlat <- c(  50,  46.1,   60,  50,  -10,  55, -10,  15)[indx]
minlon = -80.9; maxlon = -78.3; minlat = 42.4; maxlat = 44.7


# due to limited storage, let's break one year into different months
#yr  <- args[2]    # get year string from python code, YYYY e.g., '2018'
#mon <- args[3]    # get month from python code, MM, e.g., '01'

yr <- '2018'
mon <- '12'
start.date <- as.Date(paste0(yr, formatC(mon, width = 2, flag = 0), '01'), '%Y%m%d')
end.date <- as.Date(paste0(yr, formatC(mon, width = 2, flag = 0), 
                               formatC(lubridate::days_in_month(start.date),
                                       width = 2, flag = 0)), '%Y%m%d')

# timestr in form of YYYYMMDD
timestr <- gsub('-', '', seq(start.date, end.date, by = 'day'))


# ---------------------------------------------------------------------------
# paths and patterns for Tair and Tsoil files
# ---------------------------------------------------------------------------
era5.path   <- file.path(input.path, 'ERA5', yr)    # ERA5 Tair and Tsoil
daymet.path <- file.path(input.path, 'Daymetv3', yr) # for Tair
nldas.path  <- file.path(input.path, 'NLDAS', yr)    # for Tsoil

#' common portions in the filenames before YYYY* for grabbing all available files
#' here are examples of the suitable filenames: 
#' ERA5:   STL1_201801.nc (Tsoil), 2T_201801.nc (Tair)
#' Daymet: daymet_v3_tmax_2018_na.nc4, daymet_v3_tmin_2018_na.nc4
#' NLDAS:  NLDAS_NOAH0125_H.A20180101.0000.002*.nc4
#' SMUrF will search for the correct files that match @param timestr
# choose temperature products and variable names

# if you used ERA5 temp, you need to have NN models trained using era5
# nn.indx will decide which fields and NN models to use
nn.indx <- 2
TA.field   <- c('daymet', 'ERA5')[nn.indx]
TS.field   <- c('NLDAS',  'ERA5')[nn.indx]
TA.varname <- c('daymet_v3', '2T')[nn.indx]
TS.varname <- c('SoilT_0_10cm',  'STL1')[nn.indx]

# which NN model to predict Reco, these should match the temp field you chose 
# e.g., if you used ERA5 temp, you need to have NN models trained ysing era5
nn.pattern  <- c('daymet_nldas', 'era5')[nn.indx]           
nn.platform <- 'neuralnet' #'keras' 
# pretrained models are stored under "data/NN_models"
# nn.dir <- file.path(smurf_wd, 'data/NN_models/neuralnet')

# get the correct temp paths; high res daymet and NLDAS, only for US 
TA.path <- ifelse(TA.field == 'daymet', daymet.path, era5.path)
TS.path <- ifelse(TS.field == 'NLDAS', nldas.path, era5.path)


# SM import EVI data for autotrophic respiration adjustment, 03/13/2024
EVI.path <- paste0(input.path,'/MODIS_EVI/')

EVI.pattern <- 'MODIS_V061_EVI_2018_qc_extended'

# If the EVI data has not yet been calculated from the MODIS reflectance data
# this if statement will process it (including removing bad quality flag values)
if (length(grep(EVI.pattern,list.files(EVI.path)))==0){
    print("Process MODIS EVI data")
    mod_dir='C:/Users/kitty/Documents/Research/SIF/UrbanVPRM/UrbanVPRM/dataverse_files/MODIS_reflectance/MODIS_V061_GTA_AppEEARS_2018'

    mod_EVI(mod_dir,yr,minlon,maxlon,minlat,maxlat,reg.name,lc.path,lc.pattern,lc.max.yr,EVI.path,EVI.pattern)
}


# SM added code to choose reference pixel & find minimum date of reference 
#    (needs to be uncommented), 04/18/2024
# ------------------------------------------
# *** Uncomment this to choose a reference pixel over the dominant vegetation 
## type outside the urban area:
##import the land cover
#reg.ext <- raster::extent(minlon, maxlon, minlat, maxlat) # regional extent
#lc_rast <- prep.mcd12(lc.path, lc.pattern, yr, lc.max.yr, reg.name, reg.ext)
#lc_zoom <-  prep.mcd12(lc.path, lc.pattern, yr, lc.max.yr, reg.name,raster::extent(-75.55,-75.35,45.75,45.85))
##list the most common land cover types
#sort(table(values(lc_rast)),decreasing = T) 
## Choose the most common natural vegetated land cover (i.e. not urban (13), 
## cropland (12), crop-natural mosaic (14), barren (16), or water bodies (17))

## List pixels over dominant vegetation type (DBF for GTA, Montreal & Ottawa)
#dom_pxls <- which(values(lc_rast)==4)
#lc_df <- as.data.frame(rasterToPoints(lc_rast))
#lc_df[dom_pxls,] #Choose a pixel that is surrounded by that land cover type

# End of Uncomment *** ---------------------

# **** IF YOU HAVE CHANGED THE EXTENT THIS ALSO NEEDS TO BE CHANGED!!!! ****
# Select the index of a reference pixel over a pixel with the dominant 
# vegetation type outside of the urban area (deciduous forest for GTA & Montreal):
# For the GTA I used 132707
# For Montreal/Ottawa I used 129840
# For Montreal area used 17041 (-73.89792, 45.90208)
# For Ottawa area used 16044 (-75.45208, 45.80208)
EVI_index <- 132707

# ---------------------------------------------
# *** Uncomment to determine the date where EVI is a minimum at the reference pixel
#if (as.numeric(yr)%%4==0){max_day=366}else{max_day=365}
#for (i in 1:max_day){
#    if (i==1){
#        EVI_list <- raster(paste0(EVI.path,'/',EVI.pattern,'.tif'),band=i)[EVI_index]
#    }else{
#        EVI_list <- append(EVI_list,raster(paste0(EVI.path,'/',EVI.pattern,'.tif'),band=i)[EVI_index])
#    }
#}
#min_date <- which(EVI_list==min(EVI_list))
##end of determine date *** -------------------

# --------------------------
# Uncomment if you already know the date where EVI is a minimum at the reference pixel
# See code block above to determine date with minimum NEE at the reference pixel

# *** CHANGE THE MIN DATE FOR THE DAY OF THE YEAR WITH THE LOWEST EVI AT THE 
# REFERENCE PIXEL.***
min_date <- 49
# end of uncomment ---------

EVI_ref_min <- raster(paste0(EVI.path,'/',EVI.pattern,'.tif'),band=min_date)[EVI_index] 

# For 2018 in the Ottawa region (V061): 
# Lowest value at reference pixel 16044 occurs on day 86 (End of March)

# For 2018 in the Montreal region (V061): 
# Lowest value at reference pixel 17041 occurs on day 365 (December)

# For 2018 in the Ottawa/Montreal region (V061): 
  # Lowest value at reference pixel 129840 occurs on day 365 (December)

# For 2018 in the GTA (V061): 
  # Lowest value at reference pixel 132707 occurs on day 49 (February)
# For 2019 in the GTA V061:
  # Lowest value at reference pixel 132707 occurs on day 50 (February)
# For 2020 in the GTA V061:
  # Lowest value at reference pixel 132707 occurs on day 5 (January)
# For 2021 in the GTA V061: 
  # Lowest value at reference pixel 132707 occurs on day 358 (December)
# For the GTA using MODIS V006: 
  # The lowest value at the reference pixel occurs on day 5 of the year for 2018 ( I think there was a bug here)
  # and day 50 for 2019 (February), day 3 for 2020

# ---------------------------------------------------------------------------
# use SLURM for parallel simulation settings
# ---------------------------------------------------------------------------
# too many cores may slow the calculations and cause job being killed
n_nodes  <- 1
n_cores  <- 1       # max of 5 cores if running on CHPC @utah
job.time <- '24:00:00'      # total job time
slurm    <- n_nodes > 1  # logical, TF
slurm_options <- list(time = job.time, account = 'lin-kp', partition = 'lin-kp')
jobname <- paste('SMUrF_Reco', reg.name, yr, sep = '_') 
message(jobname)
#stop()

# ----------------------------------------------------------------------------
# Start running Reco model 
# ----------------------------------------------------------------------------
message('Initializing Reco estimates')
message('Number of parallel threads: ', n_nodes * n_cores)
smurf_apply(FUN = predReco_biome, slurm, slurm_options, n_nodes, n_cores, jobname, 
            reg.name, reg.path, minlon, maxlon, minlat, maxlat, timestr, 
            lc.path, lc.pattern, lc.max.yr, TA.path, TA.field, TA.varname,
            TS.path, TS.field, TS.varname, nn.pattern, nn.platform, 
            ISA.path, ISA.pattern, isa.max.yr,isa.tf,
            EVI.path, EVI.pattern, EVI_index, EVI_ref_min, smurf_wd, tmpdir)

removeTmpFiles(h=0.25) #remove temporary files older than 15 minutes

print('Done!')
q('no')

# end of script


# ----------------------------------------------------------------------------
# script to re-run SMUrF for missing time stamps, if they exist
# ----------------------------------------------------------------------------
if (F) {

    reco.path <- file.path(output.path, reg.name, 'daily_mean_Reco_neuralnet', nn.pattern, yr)
    reco.files <- list.files(reco.path, '.nc')
    exist.timestr <- substr(reco.files, nchar(reco.files) - 10, nchar(reco.files) - 3)
    
    all.timestr <- gsub('-', '', seq(as.Date(paste0(yr, '-01-01')), 
                                     as.Date(paste0(yr, '-12-31')), by = 'day'))
    miss.timestr <- all.timestr[!substr(all.timestr, 1, 8) %in% exist.timestr]
    print(miss.timestr)
    
    jobname <- paste('SMUrF_Reco', reg.name, yr, 'missing', sep = '_')
    slurm_options <- list(time = '06:00:00', account = 'lin-kp', partition = 'lin-kp')
    smurf_apply(FUN = predReco_biome, slurm = T, slurm_options, n_nodes = 8, 
                n_cores = 3, jobname, reg.name, reg.path, 
                minlon, maxlon, minlat, maxlat, timestr = miss.timestr, 
                lc.path, lc.pattern, TA.path, TA.field, TA.varname, TS.path, 
                TS.field, TS.varname, nn.pattern, nn.platform, smurf_wd, tmpdir)

    q('no')            
}   # end if

