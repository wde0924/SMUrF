# script to preprocess observed variables and interpolate modeled variables
# DW, 05/10/2019 

# extract Tsoil from NLDAS or ERA5 right at the fluxsite 
# always use the shallowest soil layer for Tsoil, 07/30/2019 

#### source all functions and load all libraries
homedir  <- '/uufs/chpc.utah.edu/common/home'
smurf_wd <- file.path(homedir, 'lin-group7/wde/SMUrF')
source(file.path(smurf_wd, 'r/dependencies.r'))

inpath  <- file.path(homedir, 'lin-group7/wde/input_data')
outpath <- file.path(homedir, 'lin-group7/wde/SMUrF_output/extract_rds/fluxnet')
nldas.path  <- file.path(inpath, 'NLDAS')       # NLDAS for Tsoil (US)
daymet.path <- file.path(inpath, 'Daymetv3')    # Daymet for Tair (US) 
era5.path   <- file.path(inpath, 'ERA5')        # ERA5 for Tsoil level 1
csif.path   <- file.path(inpath, 'CSIF/clear_sky')
#tsif.path  <- file.path(inpath, 'TROPOMI/SIF_AJT_20190902')
tsif.path   <- file.path(inpath, 'TROPOMI/SIF_AJT_20200108/data')
tmpdir <- '/scratch/local/u0947337'     # temporary path to store big rasters


# --------------------------------------------------------------------------- #
# if pre-processing FLUXNET 2015
# --------------------------------------------------------------------------- #
flux.path <- file.path(inpath, 'FLUXNET')       # path for EC tower
flux.xlsx.file <- file.path(flux.path, 'FLX_AA-Flx_BIF_LATEST.xlsx')

# load daily flux data, flux unit: gC/m2/day according to FLUXNET2015
flux.files <- list.files(flux.path, 'FULLSET_D', recursive = T, full.names = T)
flux.info  <- strsplit.to.df(basename(flux.files), sep = '_')
all.sites  <- flux.info$V2

# --------------------------------------------------------------------------- #
if (F) {  # each core will work on each site

slurm_options <- list(time = '00:20:00', account = 'lin-kp', partition = 'lin-kp')
obs.fn <- smurf_apply(FUN = extract.fluxnet.daily, slurm = T, slurm_options,
                      n_nodes = 2, n_cores = 10, jobname = 'prep_fluxnet_obs', 
                      flux.xlsx.file, site.name = all.sites, flux.path, 
                      outpath, overwriteTF = T, smurf_wd)
}   # end if


stop() 

# --------------------------------------------------------------------------- #
# *** if you finish running the previous part of the code
# run the following to extract modeled variables at flux sites
# --------------------------------------------------------------------------- #
# get rds files pre-processed by 'prep_NN_train_reco.r'
rds.files <- list.files(outpath, 'OBS_FLUXNET_', full.names = T)
rds.obs.files <- rds.files[!grepl('model', rds.files)]
obs.info <- strsplit.to.df(basename(rds.obs.files), sep = '_')
obs.sites <- obs.info$V4

# select sites that have data in and beyond 2010
sel.info  <- obs.info %>% filter(substr(V7, 1, 4) >= '2010')
sel.sites <- sel.info$V4

if (F) {

# TEMP + SIF + soil water extraction for each site
slurm_options <- list(time = '120:00:00', account = 'lin-kp', partition = 'lin-kp')
smurf_apply(FUN = extract.feature.mod, slurm = T, slurm_options, n_nodes = 4, 
            n_cores = 4, jobname = 'prep_reanalysis', outpath, 
            site.name = sel.sites, flux.source = 'FLUXNET', nldas.path, 
            daymet.path, era5.path, csif.path, csif.min.yr = 2010, tsif.path, 
            smurf_wd, tmpdir)
} 

stop() 


# --------------------------------------------------------------------------- #
# aggregate all sites' fluxes into one rds file
# --------------------------------------------------------------------------- #
# no Tsoil or not enough data
na.ts.sites <- c('US-PFa', 'BE-Bra', 'US-Wi9', 'BR-Sa1', 'FR-Fon', 'AU-Rob', 
                 'AR-Vir', 'DE-Seh', 'US-Lin', 'CN-Cng', 'CN-Du3', 'CA-Qfo', 
                 'ZA-Kru', 'ES-LgS', 'GH-Ank', 'US-Wi4', 'US-Wi0', 'IT-La2', 
                 'CA-TP2', 'US-ORv', 'CA-Obs', 'CA-Oas') 
na.swc.sites <- c('AR-Vir', 'FI-Let', 'IT-La2', 'US-GBT', 'US-Wi0', 'US-Wi4', 
                  'DE-RuR', 'RU-Sam', 'FR-Pue', 'CZ-wet', 'DE-Akm', 'DE-SfN', 
                  'DE-Spw', 'DE-Zrk', 'DK-NuF', 'DK-ZaF', 'NO-Adv', 'RU-Che', 
                  'SE-St1', 'US-Los', 'US-Myb', 'US-ORv', 'US-Tw1', 'US-Tw4', 
                  'US-WPT', 'US-Ha1', 'US-Wi3', 'ES-Ln2', 'RU-Cok', 'US-Wi6', 
                  'US-Twt', 'RU-SkP')     # sites with no SWC

rds.files <- list.files(outpath, 'daily_OBS_FLUXNET_model_', full.names = T)
site.info <- strsplit.to.df(basename(rds.files), sep = '_') 

#rds.files <- rds.files[!site.info$V1 %in% c(na.ts.sites, na.swc.sites)]
rds.files <- rds.files[!site.info$V5 %in% na.ts.sites]
site.info <- strsplit.to.df(basename(rds.files), sep = '_') 
uni.site  <- unique(site.info$V5)
uni.igbp  <- unique(site.info$V6)

func.readRDS.fluxnet <- function(x) { 
    dat <- readRDS(x); sel <- dat[, colnames(dat) %in% coln]
    print(unique(sel$site)); return(sel)
}

# choose the columns you'd like to subseta
coln <- c('site', 'igbp', 'LT', 'TIMESTAMP',  'DoY', 'CSIF_mean', 'CSIF_inst',
          'TA_F_MDS', 'TA_era5', 'TS_F_MDS_1', 'TS_era5', 'GPP', 'RECO', 'NEE')  

# convert fluxes units from gC m-2 day-1 to umol m-2 s-1 
fn.dat <- do.call(rbind, lapply(rds.files, func.readRDS.fluxnet)) %>% arrange(igbp)

ref.dat <- fn.dat %>% 
           mutate(country = substr(site, 1, 2), 
                  GPP = GPP / 12 * 1E6 / 86400, GPP_unit = 'umol m-2 s-1',
                  RECO = RECO / 12 * 1E6 / 86400, Reco_unit = 'umol m-2 s-1',
                  NEE = NEE / 12 * 1E6 / 86400, NEE_unit = 'umol m-2 s-1',
                  TA_unit = 'degC', TS_unit = 'degC') %>% rename(Reco = RECO)

write.table(ref.dat, file = file.path(smurf_wd, 'data/daily_fluxes_global.csv'), 
            sep = ',', row.names = F)

# choose the columns you'd like to subset
coln <- c('site', 'igbp', 'LT', 'TIMESTAMP',  'DoY', 'CSIF_mean', 'CSIF_inst',
          'TA_F_MDS', 'TA_era5', 'TA_daymet', 'TS_F_MDS_1', 'TS_era5', 'TS_nldas',
          'GPP', 'RECO', 'NEE')   

us.files <- rds.files[grepl('US', rds.files)]
fn.dat <- do.call(rbind, lapply(us.files, func.readRDS.fluxnet)) %>% arrange(igbp)

ref.dat <- fn.dat %>% 
           mutate(country = substr(site, 1, 2), 
                  GPP = GPP / 12 * 1E6 / 86400, GPP_unit = 'umol m-2 s-1',
                  RECO = RECO / 12 * 1E6 / 86400, Reco_unit = 'umol m-2 s-1',
                  NEE = NEE / 12 * 1E6 / 86400, NEE_unit = 'umol m-2 s-1',
                  TA_unit = 'degC', TS_unit = 'degC') %>% rename(Reco = RECO)

write.table(ref.dat, file = file.path(smurf_wd, 'data/daily_fluxes_conus.csv'), 
            sep = ',', row.names = F)

