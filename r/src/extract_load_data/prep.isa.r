#' @author Sabrina Madsen-Colford, 11/10/2022
#' This code prepares impervious surface area (ISA) data for use in SMUrF

prep.isa <- function(isa.path, isa.pattern, yr = '2018', isa.max.yr, 
                       reg.name = NULL, reg.ext) {
  
  # load all ISA data and find the correct one according to @param yr
  isa.files <- list.files(isa.path, isa.pattern, full.names = T)
  
  if (length(isa.files) > 1) {
    isa.file <- isa.files[grepl(yr, isa.files)]
    if (yr > isa.max.yr) isa.file <- isa.files[length(isa.files)]
  } else isa.file = isa.files  
  
  # make sure ISA domain > your target `reg.ext`
  if (length(isa.file) == 0) {  # if no ISAland cover data found, create shapefile
    create.shapefile(reg.name, reg.ext, outpath = isa.path, overwrite = T) 
    stop(paste('predGPP(): NO Impervious Surface Area data for required region 
               in', yr, '\nPlease check created shapefiles under', isa.path, 
               'and download local ISA data (For Toronto using the city of 
               Torontos Impermeable Surface product (https://ckan0.cf.opendata.
               inter.prod-toronto.ca/gl/dataset/topographic-mapping-impermeable
               -surface), Canadian Annual Crop Inventory (ACI) (https://www.agr.
               gc.ca/atlas/apps/metrics/index-en.html?appid=aci-iac), the 
               Southern Ontario Land Resource Information System (SOLRIS) 
               (https://www.arcgis.com/home/item.html?id=0279f65b82314121b5b5ec
               93d76bc6ba) and the Global Man-made Impervious Surface (GMIS)
               (https://gis.earthdata.nasa.gov/portal/home/item.html?id=38f233d3
               5da34e0cada61bc19faa0819)). Combine using 
               Toronto_permeability_plot_GMIS.r'))
  } else {
    
    if (grepl('nc',  isa.file)) isa.stk <- stack(isa.file, varname = 'ISA')
    if (grepl('tif', isa.file)) isa.stk <- raster(isa.file)
    
    if ( nlayers(isa.stk) > 1 ) {
      isa.names   <- as.numeric(substr(gsub('X', '', names(isa.stk)), 1, 4))
      layer.indx <- findInterval(yr, isa.names)
      isa.rt      <- raster::crop(subset(isa.stk, layer.indx), reg.ext)
    } else isa.rt <- raster::crop(isa.stk, reg.ext)   # end if subset layers
  }   # end if LC file        
  return(isa.rt)
}