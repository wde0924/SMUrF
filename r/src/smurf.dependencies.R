# little functions needed for SMUrF (SIF for Modeling Urban Fluxes)
# written by Dien Wu, 12/10/2018

# last updates, DW, 07/19/2019

# ---------------------------------------------------------------------------- #
fun.prod <- function(x, y) {return(x * y)}


# ---------------------------------------------------------------------------- #
ggdef.col <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# ---------------------------------------------------------------------------- #
writeENVI <- function(x, f, ...) {
    writeRaster(x, f, overwrite = TRUE, ...)
    cat('band names = {', paste(names(x), collapse = ','), '}', '\n', 
        file = extension(f, 'hdr'), append = TRUE)
}


# ---------------------------------------------------------------------------- #
# subroutine to split string and convert to data frame, DW, 12/20/2018
strsplit.to.df <- function(string, sep = '_', ncol = NULL, colnms = NULL) {

    library(stringr)
    ncol <- unique(str_count(string, sep) + 1)
    info <- matrix(unlist(strsplit(string, sep)), ncol = ncol, byrow = T)
    info.df <- data.frame(info, stringsAsFactors = F)

    if (is.null(colnms)) colnms <- paste0('V', seq(1, ncol, 1))
    colnames(info.df) <- colnms
    return(info.df)
}


# ---------------------------------------------------------------------------- #
# subroutine to convert df to raster, DW, 12/21/2018, referred from 
# https://stackoverflow.com/questions/19627344/how-to-create-a-raster-from-a-data-frame-in-r
df2raster <- function(df, crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') {

    library(raster)

    # create spatial points data frame
    spg <- df
    colnms <- colnames(spg)

    if ('x' %in% colnms & 'y' %in% colnms) coordinates(spg) <- ~ x + y
    if ('lon' %in% colnms & 'lat' %in% colnms) coordinates(spg) <- ~ lon + lat

    # coerce to SpatialPixelsDataFrame
    gridded(spg) <- TRUE

    # coerce to raster
    rt <- raster(spg)
    crs(rt) <- crs

    return(rt)
}

