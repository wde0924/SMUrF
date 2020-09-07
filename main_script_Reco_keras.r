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
#' 04/02/2020 use Keras API instead for NN models, no need to separate biomes
#' 04/06/2020 need to force Reco as NA for water, ice, and barren areas
#' ---------------------------------------------------------------------------

# activate virtual environment if you run the script on CHPC
#module unload python/3.5.2  #Unload the python executable from the original distribution
#source ~/VENV3.5.2/bin/activate #for bash shell

if (F) {

    # installing python cores for Keras and TensorFlow
    install.packages('keras')
    install.packages('tensorflow')
    install.packages('reticulate')
    keras::install_keras(tensorflow = '2.0.0a0', method = 'conda')
    tensorflow::install_tensorflow() 
    reticulate::py_config()    # check python configuration, need python 3

}

#args <- commandArgs(trailingOnly = TRUE)

#### source all functions and load all libraries
homedir  <- '/uufs/chpc.utah.edu/common/home'
smurf_wd <- file.path(homedir, 'lin-group7/wde/SMUrF'); setwd(smurf_wd) 
source('r/keras.dependencies.r')              # source all functions

# ---------------------------------------------------------------------------
# Paths one needs to modify 
# ---------------------------------------------------------------------------
input.path  <- file.path(homedir, 'lin-group7/wde/input_data')
output.path <- file.path(homedir, 'lin-group7/wde/output')
tmpdir      <- '/scratch/local/u0947337'    # raster operation may need temporary disk space
lc.path     <- file.path(smurf_wd, 'data/MCD12Q1')
lc.pattern  <- 'MCD12Q1.006_LC_Type1'

## NN model platform 
nn.platform <- c('keras', 'neuralnet')[1]
nn.dir      <- file.path(smurf_wd, 'data')
nn.pattern  <- 'nn_era5_mon_season_igbp.h5'
#nn.pattern <- 'NN_model_Reco_tair_tsoil_sifgpp_'


# ---------------------------------------------------------------------------
# Variables one needs to modify 
# ---------------------------------------------------------------------------
# name your region, needs to be consistent with that in main_script_GPPv2.r
indx <- 4
reg.name <- c('westernCONUS',   'easternCONUS',   'westernEurope', 'easternChina',    
              'northernAfrica', 'southernAfrica', 'NewZealand',    'LosAngeles')[indx] 
cat(paste('Working on', reg.name, '...\n'))

#' please make sure this domain is <= than the domain of MODIS land cover,
#' (minlon, maxlon, minlat, laxlat) that matches above @param reg.name
# (-125, -95,  30,  50) for western CONUS cities e.g., SLC, LA...
# ( -95, -65,  25,  50) for eastern CONUS cities for NY, Bos, Baltimore, Indy
# ( -11,  25,  35,  60) for western european cities, e.g., Madird, Paris...
# ( 105, 125,  15,  35) for southeastern China cities
# ( 160, 180, -50, -30) for New Zealand 
# (  10,  35, -35, -20) for South African cities, e.g., CapeTown, Johannesburg
#' these lat/lon should follow the order of @param reg.name
# *** too large a spatial extent may lead to memory issue, DONT DO ENTIRE GLOBE
minlon <- c(-125, -95, -11, 100,  10,  10, 160, -120)[indx]
maxlon <- c( -95, -65,  20, 125,  35,  35, 180, -115)[indx]
minlat <- c(  25,  25,  35,  20,  20, -35, -50,   33)[indx]
maxlat <- c(  50,  50,  60,  50,  35, -20, -30,   35)[indx]
#minlon = -90; maxlon = -80; minlat = 35; maxlat = 45

yr <- '2018'
time.df <- expand.grid(yr, formatC(seq(1, 12, 1), width = 2, flag = 0), 
                           formatC(seq(1, 31, 8), width = 2, flag = 0)) %>% 
           mutate(timestr = paste0(Var1, Var2, Var3))
timestr <- time.df$timestr


# ---------------------------------------------------------------------------
# settings for Tair and Tsoil data
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

# choose temperature products and their corresponding variable names
# if you used ERA5 temp, you need to have NN models trained ysing era5
# nn.indx will decide which fields and NN models to use
# daymet + NLDAS are only available for the US
nn.indx    <- 2             
TA.field   <- c('daymet',    'ERA5')[nn.indx]   
TS.field   <- c('NLDAS',     'ERA5')[nn.indx]
TA.varname <- c('daymet_v3', '2T')[nn.indx]
TS.varname <- c('TSOIL',     'STL1')[nn.indx]


# ----------------------------------------------------------------------------
# The following part requires no changes unless an error occurs 
# ----------------------------------------------------------------------------
# get the correct temp paths; high res daymet and NLDAS, only for US 
TA.path <- ifelse(TA.field == 'daymet', daymet.path, era5.path)
TS.path <- ifelse(TS.field == 'NLDAS',  nldas.path,  era5.path)

## paths for GPP, see 'main_script_GPP.r' for more details on GPP estimates
reg.path  <- file.path(output.path, reg.name)
gpp.files <- list.files(reg.path, 'daily_mean_SIF_GPP_uncert', full.names = T)
gpp.file  <- gpp.files[grepl(yr, gpp.files)]  # one GPP file per year
if (length(gpp.file) == 0) stop('NO GPP files found...see `main_script_GPP.r`\n')


# ---------------------------------------------------------------------------
# use SLURM for parallel simulation settings
# too many cores may slow the calculations and cause job being killed
# ---------------------------------------------------------------------------
n_nodes  <- 10
n_cores  <- 1               # max of 5 cores if running on CHPC @utah
job.time <- '24:00:00'      # total job time
slurm    <- n_nodes > 1     # logical, TF
slurm_options <- list(time = job.time, account = 'lin-kp', partition = 'lin-kp')
jobname  <- paste('SMUrF_Reco', reg.name, yr, sep = '_')
#stop()


# ----------------------------------------------------------------------------
# Start running Reco model 
# ----------------------------------------------------------------------------
message('Initializing Reco estimates')
message('Number of parallel threads: ', n_nodes * n_cores)

smurf_apply(FUN = predReco_biome, slurm, slurm_options, n_nodes, n_cores, jobname, 
            reg.name, reg.path, minlon, maxlon, minlat, maxlat, timestr, 
            lc.path, lc.pattern, gpp.file, TA.path, TA.field, TA.varname, 
            TS.path, TS.field, TS.varname, nn.dir, nn.pattern, nn.platform, 
            smurf_wd, tmpdir)

q('no')

# end of script
