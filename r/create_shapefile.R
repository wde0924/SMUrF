# a script to create shapefile for downloading data from App

# if no tif file found for MCD12Q1, you can generate shapefile as below
# and download data from https://lpdaacsvc.cr.usgs.gov/appeears/
# *** choose tif format, geographic projection, LC_Type1 for different years 

lc.path  <- ''
shp.path <- file.path(lc.path, 'shapefile')     # path for MCD shapefile
dir.create(shp.path, recursive = T, showWarnings = F)

# specify regional extent and create shapefiles for downloading MCD data 
minlon <- 
maxlon <- 
minlat <- 
maxlat <- 
reg.name <- 

reg.ext <- extent(minlon, maxlon, minlat, maxlat)

# create shapefile and all dependences
shp.file <- create.shapefile(site = reg.name, site.ext = reg.ext, 
                             outpath = shp.path, overwrite = T)


# use created *.zip file to download MCD12Q1 and MOD44B data in tif format 
# from https://lpdaacsvc.cr.usgs.gov/appeears/

