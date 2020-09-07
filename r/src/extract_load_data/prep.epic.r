# subroutine to load hourly SW radiation from EPIC, DW, 05/25/2020 

# need to load 4-day SW
prep.epic <- function(epic.path, timestr, site.ext, nhrs = NULL, tz = 'UTC') {

    # ------------------------------------------------------------------------ #
    # locate the files needed for `timestr` and `nhrs` (if nhrs is not NULL)
    doy.df <- get.doy.from.timestr(timestr, nday = 1, nhrs = nhrs, tz = tz)
    seq.dates <- seq(as.POSIXct(doy.df$min.timestr, tz = tz, format = '%Y%m%d%H'), 
                     as.POSIXct(doy.df$max.timestr, tz = tz, format = '%Y%m%d%H'), 
                     by = 1 * 60 * 60)    # by in second
    
    if (!is.null(nhrs)) seq.dates <- seq.dates[-length(seq.dates)]
    seq.timestr <- sort(gsub('-', '', unique(substr(seq.dates, 1, 10))))

    # select files
    epic.files <- list.files(epic.path, '.nc', full.names = T)
    epic.timestr <- (strsplit.to.df(gsub('.nc', '', basename(epic.files))))$V5
    epic.file <- epic.files[epic.timestr %in% seq.timestr]

    # lat, lon are flipped, unit of PAR is W m-2
    epic.ext <- extent(10 * (90 - site.ext[4]), 10 * (90 - site.ext[3]), 
                       10 * (180 + site.ext[1]), 10 * (180 + site.ext[2]))
    
    epic.sw <- NULL
    for (e in 1 : length(epic.file)) {

        print(seq.timestr[e])
        sw.dir.stk <- crop(stack(epic.file[e], varname = 'sw_direct'), epic.ext)
        sw.dif.stk <- crop(stack(epic.file[e], varname = 'sw_diffuse'), epic.ext)
        sw.tot.stk <- sw.dir.stk + sw.dif.stk 

        # hours since 2015-01-01 00:30:00, 
        # initially, dates is centered in an hour interval
        # now move to start of the hour interval
        epic.dates <- as.POSIXct(as.numeric(gsub('X', '', names(sw.dir.stk))) * 3600, 
                                 tz = 'UTC', origin = '2015-01-01 00:30:00') - 30 * 60

        # need to flip raster by latitude and then transpose
        sw.tot.tf <- flip(t(sw.tot.stk), direction = 'x')
        extent(sw.tot.tf) <- site.ext
        crs(sw.tot.tf) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
        names(sw.tot.tf) <- epic.dates

        if (F) {
            # sanity check for plotting 
            s1 <- levelplot(sw.tot.tf, maxpixel = 2e6, layout = c(4, 6))
            ggsave(arrangeGrob(s1), filename = paste0('../EPIC_SW_', timestr, '.png'), 
                   width = 8, height = 12)
        }
      
        if (e == 1) epic.sw <- sw.tot.tf#; epic.par <- par.tot.tf } 
        if (e >  1) epic.sw <- stack(epic.sw, sw.tot.tf) #;epic.par <- stack(epic.par, par.tot.tf) }
    }

    return(epic.sw = epic.sw)

}



# --------------------------------------------------------------------------- #
# also grab PAR from EPIC
prep.epic.all <- function(epic.path, timestr, site.ext, nhrs = NULL, tz = 'UTC') {

    # locate the files needed for `timestr` and `nhrs` (if nhrs is not NULL)
    doy.df <- get.doy.from.timestr(timestr, nday = 1, nhrs = nhrs, tz = tz)
    seq.dates <- seq(as.POSIXct(doy.df$min.timestr, tz = tz, format = '%Y%m%d%H'), 
                     as.POSIXct(doy.df$max.timestr, tz = tz, format = '%Y%m%d%H'), 
                     by = 1 * 60 * 60)    # by in second
    
    if (!is.null(nhrs)) seq.dates <- seq.dates[-length(seq.dates)]
    seq.timestr <- sort(gsub('-', '', unique(substr(seq.dates, 1, 10))))

    # select files
    epic.files <- list.files(epic.path, '.nc', full.names = T)
    epic.timestr <- (strsplit.to.df(gsub('.nc', '', basename(epic.files))))$V5
    epic.file <- epic.files[epic.timestr %in% seq.timestr]

    # lat, lon are flipped, unit of PAR is W m-2
    epic.ext <- extent(10 * (90 - site.ext[4]), 10 * (90 - site.ext[3]), 
                       10 * (180 + site.ext[1]), 10 * (180 + site.ext[2]))
    
    epic.sw <- NULL
    for (e in 1 : length(epic.file)) {

        print(seq.timestr[e])
        sw.dir.stk <- crop(stack(epic.file[e], varname = 'sw_direct'), epic.ext)
        sw.dif.stk <- crop(stack(epic.file[e], varname = 'sw_diffuse'), epic.ext)
        sw.tot.stk <- sw.dir.stk + sw.dif.stk 

        # hours since 2015-01-01 00:30:00, 
        # initially, dates is centered in an hour interval
        # now move to start of the hour interval
        epic.dates <- as.POSIXct(as.numeric(gsub('X', '', names(sw.dir.stk))) * 3600, 
                                 tz = 'UTC', origin = '2015-01-01 00:30:00') - 30 * 60

        # need to flip raster by latitude and then transpose
        sw.tot.tf <- flip(t(sw.tot.stk), direction = 'x')
        extent(sw.tot.tf) <- site.ext
        crs(sw.tot.tf) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
        names(sw.tot.tf) <- epic.dates

        par.dir.stk <- crop(stack(epic.file[e], varname = 'par_direct'), epic.ext) 
        par.dif.stk <- crop(stack(epic.file[e], varname = 'par_diffuse'), epic.ext)
        par.tot.stk <- par.dir.stk + par.dif.stk 
        par.tot.tf  <- flip(t(par.tot.stk), direction = 'x') 
        extent(par.tot.tf) <- site.ext 
        crs(par.tot.tf) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
        names(par.tot.tf) <- epic.dates

        if (e == 1) {epic.sw <- sw.tot.tf; epic.par <- par.tot.tf } 
        if (e >  1) { epic.sw <- stack(epic.sw, sw.tot.tf) ;epic.par <- stack(epic.par, par.tot.tf) }
    }

    return(list(epic.sw = epic.sw, epic.par = epic.par))

}

