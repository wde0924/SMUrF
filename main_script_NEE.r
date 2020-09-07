#' Main script to generate net ecosystem exchange
#' @author: Dien Wu, 07/12/2019
#' ---------------------------------------------------------------------------

#' @GeneralIdeas:
#' 1. Need GPP and Reco generated by 'main_script_GPP.r' and 'main_script_Reco.r'
#' 2. Need ERA5' surface insolation and Tair for downscaling GPP and Reco

#' @Updates by Dien Wu: 
#' 10/28/2019 incorporate tree and non-tree vegetation fraction from MOD44B
#' ---------------------------------------------------------------------------

#args <- commandArgs(trailingOnly = TRUE)

#### source all functions and load all libraries
homedir <- '/uufs/chpc.utah.edu/common/home'
smurf_wd <- file.path(homedir, 'lin-group7/wde/SMUrF'); setwd(smurf_wd) 
source('r/dependencies.r')              # source all functions

# ---------------------------------------------------------------------------
# Paths one needs to modify 
# ---------------------------------------------------------------------------
# input and output paths
input.path  <- file.path(homedir, 'lin-group7/wde/input_data')
output.path <- file.path(homedir, 'lin-group7/wde/output')


# ---------------------------------------------------------------------------
# Variables one needs to modify 
# ---------------------------------------------------------------------------
# name your region, needs to be consistent with that in main_script_GPPv2.r
indx <- 2
#indx <- as.numeric(args[1])
reg.name <- c('westernCONUS', 'easternCONUS',     'westernEurope', 
              'easternChina', 'easternAustralia', 'easternAsia', 
              'southAmerica', 'centralAfrica')[indx]   
reg.path <- file.path(output.path, reg.name)

# the directory that stores daily mean Reco nc files
reco.dir <- 'daily_mean_Reco_neuralnet/era5'

# please make sure this domain is <= than the domain of MODIS land cover,
# 'minlon maxlon, minlat, laxlat' should matche the order of 'reg.name' above
# *** too large a spatial extent may lead to memory issue, DONT DO ENTIRE GLOBE
minlon <- c(-125, -95,  -11, 100,  130, 125, -65, -10)[indx]
maxlon <- c( -95, -65,   20, 125,  155, 150, -40,  20)[indx]
minlat <- c(  25,  25,   35,  20,  -40,  30, -40, -10)[indx]
maxlat <- c(  50,  50,   60,  50,  -10,  55, -10,  15)[indx]

# each processor works on each month
yr <- 2018
mons <- seq(1, 12)
#yr <- as.numeric(args[2])
#mons <- as.numeric(args[3])


## paths and variable names for loading hourly Tair and SW rad on surface data
TA.field   <- 'ERA5'
SSRD.field <- c('ERA5', 'EPIC')[1]
TA.path    <- file.path(input.path, TA.field, yr) 
SSRD.path  <- file.path(input.path, SSRD.field, yr) 

# common portions in the filenames before YYYY* for grabbing all available files
# here are examples of the suitable filenames: 
# ERA5: STL1_201801.nc (Tsoil), 2T_201801.nc (Tair), SSRD_201801.nc
#' the model itself will find the correct files that match @param timestr
# no need of variable names if using EPIC--
# SMUrF will grab both direct and diffuse SW rad
TA.varname   <- '2T'
SSRD.varname <- c('SSRD', NA)[1]    


# ----------------------------------------------------------------------------
# parallel calculations
# ----------------------------------------------------------------------------
# use SLURM for parallel simulation settings
jobtime <- '1:00:00'      # total job time
n_nodes <- 1   
n_cores <- 1
slurm   <- T      # False for not running parallelly
slurm_options <- list(time = jobtime, account = 'lin-kp', partition = 'lin-kp')
jobname <- paste('SMUrF_NEE', reg.name, yr, sep = '_') 
message(jobname)
#stop()

# ----------------------------------------------------------------------------
# Start running and storing hourly fluxes in nc files (by months)
# ----------------------------------------------------------------------------
message('\n\nInitializing hourly NEE estimates')
message('Number of parallel threads: ', n_nodes * n_cores)
all.yyyymm <- paste0(yr, formatC(mons, width = 2, flag = 0))
smurf_apply(FUN = predNEE, slurm, slurm_options, n_nodes, n_cores, jobname, 
            reg.name, reg.path, reco.dir, yyyymm = all.yyyymm, TA.path, TA.field, 
            TA.varname, SSRD.path, SSRD.field, SSRD.varname, smurf_wd)

q('no')
# end of script