#' Main script to generate downscaled fluxes using 250m fveg from MODIS VCF
#' only perform this for nearfield area (i.e., region around target city)
#' @author: Dien Wu, 11/30/2019
#' ---------------------------------------------------------------------------

#### source all functions and load all libraries
homedir <- '/uufs/chpc.utah.edu/common/home'
smurf_wd <- file.path(homedir, 'lin-group7/wde/SMUrF'); setwd(smurf_wd) 
source('r/dependencies.r')              # source all functions
register_google('')     # insert your google API

#' ---------------------------------------------------------------------------
#' Variables one needs to modify 
#' ---------------------------------------------------------------------------
# name your region, needs to be consistent with that in main_script_GPPv2.r
indx <- 8
site <- c('Guangzhou', 'Shanghai',   'Beijing',       'Wuhan',        'Taipei', 
          'NewYork',   'Chicago',    'Indianapolis',  'Boston',       'Baltimore', 
          'Dallas',    'LosAngeles', 'SaltLakeCity',  'SanFrancisco', 'Phoenix', 
          'London',    'Paris',      'Berlin',        'Madrid',       'Milan')[indx]
reg.name <- c(rep('easternChina', 5), rep('easternCONUS', 5), 
              rep('westernCONUS', 5), rep('westernEurope', 5))[indx]  

# each processor works on a day in a year 
min.datestr <- '2017-06-01'
max.datestr <- '2017-09-30'

# a box centered around your city, dlat = dlon = 2 meaning 4x4degree box 
dlat <- 2.5     # in degrees
dlon <- 2.5  
loc  <- get.lon.lat(site, dlon, dlat)    # return city lon/lat and other info
# e.g., > loc
#          site  sitelon  sitelat             tz                countryid
#1 SaltLakeCity -111.891 40.76078 America/Denver United States of America
#          regid iso3   minlon   maxlon   minlat   maxlat
#1 North America  USA -114.391 -109.391 38.26078 43.26078

# get min/max lat/lon for downscaling fluxes, must be rounded to an integer
minlon <- floor(loc$minlon); maxlon <- ceiling(loc$maxlon)
minlat <- floor(loc$minlat); maxlat <- ceiling(loc$maxlat)


#' ---------------------------------------------------------------------------
#' Paths one needs to modify 
#' ---------------------------------------------------------------------------
#' raster operations may need temporary disk space for large calculations
#tmpdir <- '/scratch/local/u0947337'
tmpdir <- NA 

# input and output paths
input.path  <- file.path(homedir, 'lin-group7/wde/input_data')
output.path <- file.path(homedir, 'lin-group7/wde/output')
reg.path    <- file.path(output.path, reg.name)   # NEE path


# path for MOD44B VCF data, tree and non-tree veg fractions
# pre-processed from https://lpdaacsvc.cr.usgs.gov/appeears/
# required files: MOD44B.006.Percent_NonVegetated_doy*_<vcf.pattern>
#                 MOD44B.006.Percent_Tree_Cover_doy*_<vcf.pattern>
vcf.path <- file.path(input.path, 'MOD44B')
vcf.pattern <- 'CONUS.tif'


# ----------------------------------------------------------------------------
# parallel calculations
# ----------------------------------------------------------------------------
#' create daily string based on min/max.datestr for parallel computing 
dates <- seq(as.Date(min.datestr), as.Date(max.datestr), by = 'day')
all.timestr <- gsub('-', '', dates)

# use SLURM for parallel simulation settings
# each thred will work on one timestr
jobtime <- '01:00:00'      # total job time
n_nodes <- 8   
n_cores <- 6
slurm   <- T      # False for not running parallelly
slurm_options <- list(time = jobtime, account = 'lin-kp', partition = 'lin-kp')
jobname <- paste0('SMUrF_ds_', site) 
#stop()


# ----------------------------------------------------------------------------
# Start running and storing hourly fluxes in nc files (by months)
# ----------------------------------------------------------------------------
message('\n\nInitializing hourly NEE estimates')
message('Number of parallel threads: ', n_nodes * n_cores)
smurf_apply(FUN = predNEE.downscale, slurm, slurm_options, n_nodes, n_cores, jobname, 
            site, reg.path, timestr = all.timestr, minlon, maxlon, minlat, maxlat, 
            vcf.path, vcf.pattern, smurf_wd, tmpdir)

q('no')
# end of script



