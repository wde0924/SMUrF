#' script to preprocess hourly soil temp and moisture from NLDAS, DW, 05/23/2019 
#' return hourly soil temp or moisture as rasterBrick (if multiple hours) 
#'                                     or rasterLayer (if one single hour)
#' no daily mean calculation performed, return hourly variables

# NLDAS data is downloaded from subset of initial data from 
# https://disc.gsfc.nasa.gov/datasets/NLDAS_NOAH0125_H_V002/summary?keywords=NLDAS

# for this analysis, only nldas and soil moisture have been downloaded 
#' @input:
#' @param nldas.pattern the common portion of all the NLDAS files, before YYYYMM
#' @param timestr in the form of YYYYMMDDHH 
#' @param meanTF whether to compress hourly nldas to one mean rasterLayer
#' @param nldas.level  level 1 for 0-10cm, level 2 for 10-40cm
#'                     level 3 for 40cm - 1m, level 4 for 1-2m 

prep.soil.nldas <- function(nldas.path, nldas.varname = c('TSOIL', 'SOILM')[1], 
                            timestr = '2015061820', ext = NULL, nhrs = NULL, 
                            soil.level = c(1, 2, 3, 4)[1], 
                            nldas.pattern = 'NLDAS_NOAH0125_H.A', tz = 'UTC') {
            
    library(raster); library(dplyr)

    # list all NLDAS files and get its info
    # example file name: 'NLDAS_NOAH0125_H.A20181231.1900.002.*.nc4'
    nldas.files <- list.files(nldas.path, nldas.pattern, full.names = T, recursive = T)
    #print(str(nldas.files))
    if (length(nldas.files) == 0) stop('prep.soil.nldas(): NO NLDAS file found...\n')
    
    nldas.str  <- str_remove(basename(nldas.files), nldas.pattern)
    nldas.info <- strsplit.to.df(nldas.str, sep = '\\.') %>% 
                  mutate(timestr = paste0(substr(V1, 1, 8), substr(V2, 1, 2))) 

    # locate the files needed for `timestr` and `nhrs` (if nhrs is not NULL)
    nldas.doy.df <- get.doy.from.timestr(timestr, nday = 1, nhrs = nhrs, tz = tz)
    nldas.indx <- which(nldas.info$timestr <= nldas.doy.df$max.timestr & 
                        nldas.info$timestr >= nldas.doy.df$min.timestr)
    if (length(nldas.indx) == 0) stop(paste('NO NLDAS file found for', timestr))
    
    # select qualified files and read all in as rasterBrick
    nldas.file     <- nldas.files[nldas.indx]
    sel.nldas.info <- nldas.info[nldas.indx, ]

    # read all required hourly nldas, convert deg K to deg C
    cat(paste('prep.soil.nldas(): Take a while to read hourly', nldas.varname, 
              'from multiple files...\n'))
    # level 1 for 0-10cm, level 2 for 10-40cm
    # level 3 for 40cm - 1m, level 4 for 1-2m 
    nldas.brk <- do.call(stack, lapply(nldas.file, function(file) 
                                       raster(file, varname = nldas.varname, 
                                                    level = soil.level)))                                      

    # and crop the data according to ext (can be generated from extent())
    if (!is.null(ext)) nldas.brk <- crop(nldas.brk, ext)
    names(nldas.brk) <- sel.nldas.info$timestr  # assign names to layers

    # unit corrections
    if (nldas.varname == 'TSOIL') {
        cat('prep.soil.nldas(): *** unit corrections (degC for TSOIL)\n')
        nldas.brk <- calc(nldas.brk, fun = function(x) { x - 273.15 }) 
    }

    return(nldas.brk)
}   # end of subroutine