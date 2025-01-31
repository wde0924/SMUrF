# subroutine to grab 0.5deg GPP measurements from FLUXCOM
# Dien Wu, 04/27/2019

# `ext`: can be generated by extent() using raster package, 
#        if NULL, no need to crop the spatial domain of gridded GPP
# `timestr` and `nhrs`: timestr of YYYYMMDDHH indicating desired dates, 
#                       nhrs of numbers of hours away from timestr, can be +/-
#                       if both are not given, no need to select days

grab.fluxcom <- function(fc.path, fc.file, ext = NULL, timestr = NULL, 
                         nhrs = NULL) {

    rt <- stack(file.path(fc.path, fc.file))
    if (!is.null(ext)) rt <- crop(rt, ext)

    if (!is.null(timestr)) {

        # find the correct doy based on timestr
        doy.info <- get.doy.from.timestr(timestr, nday = 1, nhrs)
        doy.vec  <- seq(min(doy.info$find.doy, doy.info$find.doy2), 
                        max(doy.info$find.doy, doy.info$find.doy2))

        rt <- subset(rt, doy.vec)
    }

    return(rt)
}