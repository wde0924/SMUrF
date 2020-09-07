#' subroutine to calculate GPP and its uncertainty from SIF
#' @author Dien Wu, 06/10/2019

#' @param sif.rt rasterLayer, 0.05deg spatial SIF
#' @param prep.stk rasterStack, names of "SLP.AGB"  "SLP.mean" "SLP.sd"
#' 04/07/2020, DW, no need to perform unit conversion to GPP 

compute.gpp <- function(sif.rt, prep.stk) {

    # calculate the best estimated GPP (mean values) and GPP-SIF slopes 
    # slopes are gap filled for urban using AGB and weighted mean approaches
    if (raster::compareRaster(sif.rt, prep.stk) == F) 
        stop('SIF and slopes grids have different extent...
              one possibility: a larger target extent than MCD12Q1"s extent...
              please check\n')
    
    # SLOPE * SIF + INT = GPP in umol m-2 s-1
    gpp.mean <- raster::overlay(prep.stk$SLP, sif.rt, fun = fun.prod) + prep.stk$INT
    gpp.sd  <- gpp.mean * prep.stk$GPP_CV
    gpp.stk <- stack(gpp.mean, gpp.sd); names(gpp.stk) <- c('GPP_mean', 'GPP_sd')
    gpp.stk[gpp.stk < 0] <- 0 
    
    return(gpp.stk)
}   # end of subroutine
