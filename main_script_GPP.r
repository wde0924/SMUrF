#' Main script to generate daily mean SIF-based GPP for each yr
#' grab spatial SIF, assign GPP-SIF slopes with gap filling for urban core

#' @author: Dien Wu, 04/19/2019
#' last update, 03/28/2020
#' ---------------------------------------------------------------------------

#' @datasets required include: 
#' 1. OCO-2 SIF and 0.05 degree CSIF
#' 2. MODIS land cover (500m MCD12) downloaded from AρρEEARS which reforms MODIS 
#'    product to desired format,  https://lpdaacsvc.cr.usgs.gov/appeears/
#'    This link requires shapefiles for target region, which can be generated 
#'    using `create.shapefile.r`
#' 3. Aboveground biomass (GlobBiomass) downloaded from 
#'    http://globbiomass.org/wp-content/uploads/GB_Maps/Globbiomass_global_dataset.html
#' ---------------------------------------------------------------------------

#' @updates by DW:
#' 04/19/2019: rearrange and optimize all subroutines. 
#' 04/22/2019: generate temporal GPP fluxes
#' 05/21/2019: optimize the code
#' 05/30/2019: store as RasterStack instead of data.frame in .envi files
#' 06/10/2019: add uncertainty to GPP-SIF slopes from Zhang et al., 2018
#'             add an alternative way to gap fill urban slopes (mean of slopes)
#' 06/12/2019: speed up the calculation, by getting the GPP-SIF slopes first, 
#'             before looping over each day; also divide CONUS to east vs. west
#' 07/03/2019, get rid of the AGB-scaled GPP-SIF slopes, only use the weighted mean slopes
#' 07/03/2019, add parallel computation to GPP estimates, one core -> one yr
#' 07/25/2019, use GOSIF (every 8 days) and CSIF (every 4 days) as errors for SIF
#' 07/26/2019, store all daily GPP into one nc file, see save.raster2nc.R
#' 09/05/2019, gap filling for all cities within a spatial extent
#' 03/29/2020, update urban gap-fill and incorporate C3-C4 partitioning
#' 04/07/2020, use GPP-CSIF slopes based on GPP in units of umol m-2 s-1
#' ---------------------------------------------------------------------------

# source all functions and load all libraries
homedir <- '/uufs/chpc.utah.edu/common/home'
smurf_wd <- file.path(homedir, 'lin-group7/wde/SMUrF'); setwd(smurf_wd)
source('r/dependencies.r') 

# ---------------------------------------------------------------------------
# Paths one needs to modify 
# ---------------------------------------------------------------------------
# input: e.g., OCO-2, spatial SIF, above ground biomass
input.path  <- file.path(homedir, 'lin-group7/wde/input_data')
output.path <- file.path(homedir, 'lin-group7/wde/output')

# path for spatial CSIF, Zhang et al., 2018
csif.cpath <- file.path(input.path, 'CSIF/clear_sky')   # clearsky CSIF

# path for 100m AGB from GlobBiomass, need to download 40x40deg tiles of data
# from http://globbiomass.org/wp-content/uploads/GB_Maps/Globbiomass_global_dataset.html
agb.path <- file.path(input.path, 'biomass/GlobBiomass')

# path for 500m IGBP, need to download from https://lpdaacsvc.cr.usgs.gov/appeears/
lc.path    <- file.path(smurf_wd, 'data/MCD12Q1')
lc.pattern <- 'MCD12Q1.006_LC_Type1'

# indicate the latest year available of MCD12Q1
# if no data beyond 2018, use 2018 LC for 2019 and beyond
lc.max.yr <- 2018  
lc.res    <- 1/240     # horizontal grid spacing of land cover in degrees

# raster operations may need temporary disk space for large calculations
#tmpdir <- '/scratch/local/u0947337'
tmpdir <- NA

# ---------------------------------------------------------------------------
# Variables one needs to modify 
# ---------------------------------------------------------------------------
# name and choose your desired region/nation and set a spatial extent
indx <- 3
reg.name <- c('westernCONUS', 'easternCONUS',     'westernEurope', 
              'easternChina', 'easternAustralia', 'easternAsia', 
              'southAmerica', 'centralAfrica')[indx]   
cat(paste('Working on', reg.name, '...\n'))

# output paths to store GPP results
gpp.path <- file.path(output.path, reg.name)
dir.create(gpp.path, recursive = T, showWarnings = F)

# make sure the domain to be modeled is <= than the domain of MODIS land cover
#' (minlon, maxlon, minlat, laxlat) that matches above @param reg.name
#' these lat/lon should follow the order of @param reg.name
# *** too large a spatial extent may lead to memory issue, DONT DO ENTIRE GLOBE
minlon <- c(-125, -95,  -11, 100,  130, 125, -65, -10)[indx]
maxlon <- c( -95, -65,   20, 125,  155, 150, -40,  20)[indx]
minlat <- c(  25,  25,   35,  20,  -40,  30, -40, -10)[indx]
maxlat <- c(  50,  50,   60,  50,  -10,  55, -10,  15)[indx]
#minlon = -90; maxlon = -80; minlat = 35; maxlat = 45

# *** choose yrs, if multiple years, each thred will work on one year
all.yrs <- seq(2010, 2014)

# ----------------------------------------------------------------------------
# specify CSIF related parameters
# ----------------------------------------------------------------------------
sif.prod <- 'CSIFclear'
sif.var  <- 'clear_daily_SIF'
sif.nd   <- 4         # temporal resoultion, every 4 or 8 days
sif.res  <- 0.05      # 0.05 deg res
sif.rmTF <- TRUE      # if TRUE, force negative SIF values as zero
sif.path <- csif.cpath 

# txtfile that store GPP-SIF relation
slp.file <- file.path(smurf_wd, 'data/GPP-CSIF_stat_Wu.txt') 

#' whether to re-generate and store 500 m GPP-SIF slopes as tif files in 'gpp.path'
#' if you change @param slp.file, you need to turn on this flag
slp.overwriteTF <- T


# ---------------------------------------------------------------------------
# slurm settings, each core will work on one year of GPP
# ---------------------------------------------------------------------------
n_nodes <- 5
n_cores <- 1                # max core # of 2 to avoid jobs being killed
job.time <- '24:00:00'      # total job time
slurm <- T                  # logical, whether to run parallel-ly
slurm_options <- list(time = job.time, account = 'lin-kp', partition = 'lin-kp')
jobname <- paste0('SMUrF_GPP_', reg.name)

# ---------------------------------------------------------------------------
# read and process above ground biomass data before parallel computing
# ---------------------------------------------------------------------------
cat('Preparing 100 m AGB and 500 m C3-C4 ratio before submitting jobs...\n')
agb.file <- grab.agb.GlobBiomass(agb.path, minlon, maxlon, minlat, maxlat) 

# load 10 km x 10 km C3-C4 ratio using pre-processed nc file, DW, 03/29/2020
# project 10 km ratios to 500 m that matches MCD12
ratio.file <- prep.c4.ratio(smurf_wd, lc.path, lc.pattern, yr = all.yrs[1], 
                            lc.max.yr, reg.name, minlon, maxlon, minlat, maxlat, 
                            gpp.path)

# ----------------------------------------------------------------------------
# Start running SIF-based GPP model 
# ----------------------------------------------------------------------------
message('Initializing GPP estimates')
message('Number of parallel threads: ', n_nodes * n_cores)

# start running model parallel-ly, only `all.yrs` can be a vector
# other variables should be a vector or a list
smurf_apply(FUN = predGPP, slurm, slurm_options, n_nodes, n_cores, jobname, 
            reg.name, minlon, maxlon, minlat, maxlat, yr = all.yrs, 
            lc.path, lc.pattern, lc.max.yr, lc.res, agb.file, ratio.file, 
            sif.prod, sif.path, sif.var, sif.nd, sif.res, sif.rmTF, 
            gpp.path, smurf_wd, slp.file, slp.overwriteTF, tmpdir)

q('no')

### end of script








# ----------------------------------------------------------------------------
# short script for plotting gap-fileed Land cover, SIF, and GPP 
# ----------------------------------------------------------------------------
homedir <- '/uufs/chpc.utah.edu/common/home'
smurf_wd <- file.path(homedir, 'lin-group7/wde/SMUrF'); setwd(smurf_wd)
source('r/dependencies.r'); library(gridExtra); library(rasterVis)
#register_google(key = '')

# choose a day and cities 
timestr <- '20160617'
reg <- c('westernCONUS', 'easternCONUS', 'Europe')[2]
reg.ext <- c('-125_-95_30_50', '-95_-65_25_50')[2]
if (reg == 'westernCONUS') sites <- c('LosAngeles', 'SaltLakeCity', 'Seattle')
if (reg == 'easternCONUS') sites <- c('Baltimore', 'Boston', 'Indianapolis')
yr <- substr(timestr, 1, 4)

# GPP and gap-filled MCD path and file
gpp.path <- file.path(homedir, 'lin-group7/wde/output', reg)
gpp.files <- list.files(gpp.path, 'daily_mean_SIF_GPP_uncert', full.names = T)
lc.path  <- file.path(homedir, 'lin-group7/wde/input_data/MCD12Q1')
lc.files  <- list.files(lc.path, reg.ext, full.names = T)

# read all data
gpp.file <- gpp.files[grepl(yr, gpp.files)]
lc.file  <- lc.files[grepl(yr, lc.files)]
gpp.stk  <- stack(gpp.file, varname = 'GPP_mean')
sif.stk  <- stack(gpp.file, varname = 'SIF_mean')
gpp.sec  <- as.numeric(gsub('X', '', names(gpp.stk)))
all.dates <- as.POSIXct(gpp.sec, origin = '1970-01-01 00:00:00', tz = 'UTC')
all.timestr <- format(all.dates, '%Y%m%d')	# YYYYMMDD
gpp.rt <- subset(gpp.stk, which(all.timestr == timestr))   # choose a day 
sif.rt <- subset(sif.stk, which(all.timestr == timestr))   # choose a day 
lc.rt  <- raster(lc.file)

### plotting 
source('r/dependencies.r') 
plotlist <- list()
col <- rasterTheme(region = rev(terrain.colors(12)))
for (s in 1: length(sites)) {
    loc <- suppressWarnings(get.lon.lat(sites[s], dlat = 2, dlon = 2))
    ext <- extent(loc$minlon, loc$maxlon, loc$minlat, loc$maxlat)
    tmp.sif <- crop(sif.rt, ext); names(tmp.sif) <- 'CSIF'
    tmp.lc <- crop(lc.rt, ext); tmp.gpp <- crop(gpp.rt, ext)

    map <- ggplot.map(map = 'ggmap', color = 'bw', zoom = 8, 
                      center.lat = loc$sitelat, center.lon = loc$sitelon)
    s1 <- ggmap.csif(tmp.sif, sites[s], timestr, map)
    g1 <- ggmap.gpp(tmp.gpp, sites[s], timestr, map)
    l1 <- ggmap.igbp(tmp.lc, sites[s], yr, map, res = 'Gap-filled 500 m')
    lsg <- ggarrange(l1, s1, g1, nrow = 3, heights = c(1, 1, 0.95))
    plotlist[[s]] <- lsg
}   # end for s

merge.plot <- ggarrange(plotlist = plotlist, ncol = length(sites), nrow = 1) 
fn <- file.path(homedir, 'lin-group7/wde/paper3', 
                paste0('LC_SIF_GPP_', reg, '_zoomed_', timestr, '.png'))
ggsave(merge.plot, filename = fn, width = 12, height = 15)

