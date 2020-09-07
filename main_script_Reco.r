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

# when using runthem.py, turn this on
#args <- commandArgs(trailingOnly = TRUE)

#### source all functions and load all libraries
homedir <- '/uufs/chpc.utah.edu/common/home'
smurf_wd <- file.path(homedir, 'lin-group7/wde/SMUrF'); setwd(smurf_wd) 
source('r/dependencies.r')              # source all functions


# ---------------------------------------------------------------------------
# Paths one needs to modify 
# ---------------------------------------------------------------------------
input.path  <- file.path(homedir, 'lin-group7/wde/input_data')
output.path <- file.path(homedir, 'lin-group7/wde/output')

## path for the updated 500m IGBP generated from main_script_GPP.r
lc.path <- file.path(smurf_wd, 'data/MCD12Q1')
lc.pattern <- 'MCD12Q1.006_LC_Type1'
#tmpdir <- '/scratch/local/u0947337'  # raster operation may need temporary disk space
tmpdir <- NA

# ---------------------------------------------------------------------------
# Variables one needs to modify 
# ---------------------------------------------------------------------------
# name your region, needs to be consistent with that in main_script_GPPv2.r
#indx <- as.numeric(args[1])     # get region indx from python code, e.g., 1
indx <- 1
reg.name <- c('westernCONUS', 'easternCONUS',     'westernEurope', 
              'easternChina', 'easternAustralia', 'easternAsia', 
              'southAmerica', 'centralAfrica')[indx]   
              
# output path for the target region
reg.path <- file.path(output.path, reg.name)
cat(paste('Working on', reg.name, '...\n'))

# please make sure this domain is <= than the domain of MODIS land cover,
# 'minlon maxlon, minlat, laxlat' should matche the order of 'reg.name' above
# *** too large a spatial extent may lead to memory issue, DONT DO ENTIRE GLOBE
minlon <- c(-125, -95,  -11, 100,  130, 125, -65, -10)[indx]
maxlon <- c( -95, -65,   20, 125,  155, 150, -40,  20)[indx]
minlat <- c(  25,  25,   35,  20,  -40,  30, -40, -10)[indx]
maxlat <- c(  50,  50,   60,  50,  -10,  55, -10,  15)[indx]

# due to limited storage, let's break one year into different months
#yr  <- args[2]    # get year string from python code, YYYY e.g., '2018'
#mon <- args[3]    # get month from python code, MM, e.g., '01'

yr <- '2010'
mon <- '01'
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

# if you used ERA5 temp, you need to have NN models trained ysing era5
# nn.indx will decide which fields and NN models to use
nn.indx <- 2
TA.field   <- c('daymet', 'ERA5')[nn.indx]
TS.field   <- c('NLDAS',  'ERA5')[nn.indx]
TA.varname <- c('daymet_v3', '2T')[nn.indx]
TS.varname <- c('TSOIL',  'STL1')[nn.indx]

# which NN model to predict Reco, these should match the temp field you chose 
# e.g., if you used ERA5 temp, you need to have NN models trained ysing era5
nn.pattern  <- c('daymet_nldas', 'era5')[nn.indx]           
nn.platform <- 'neuralnet'
# pretrained models are stored under "data/NN_models"
# nn.dir <- file.path(smurf_wd, 'data/NN_models/neuralnet')

# get the correct temp paths; high res daymet and NLDAS, only for US 
TA.path <- ifelse(TA.field == 'daymet', daymet.path, era5.path)
TS.path <- ifelse(TS.field == 'NLDAS', nldas.path, era5.path)


# ---------------------------------------------------------------------------
# use SLURM for parallel simulation settings
# ---------------------------------------------------------------------------
# too many cores may slow the calculations and cause job being killed
n_nodes  <- 7
n_cores  <- 3       # max of 5 cores if running on CHPC @utah
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
            lc.path, lc.pattern, TA.path, TA.field, TA.varname, TS.path, 
            TS.field, TS.varname, nn.pattern, nn.platform, smurf_wd, tmpdir)

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

