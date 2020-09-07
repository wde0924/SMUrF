# compare two spatial grids 

calc.spaef <- function(x1, x2) {
    library(pracma)

    # CORR
    cc <- cor(x1, x2)

    # coefficient of variation
    cv_sim <- sd(x1, na.rm = T) / mean(x1, na.rm = T)    # CSIF
    cv_obs <- sd(x2, na.rm = T) / mean(x2, na.rm = T)    # TSIF
    alpha <- cv_sim/cv_obs

    # HISTO match
    norm.x1 <- (x1 - mean(x1))/sd(x1); norm.x2 <- (x2 - mean(x2))/sd(x2)

    # False makes Density instead of frequency, try it
    h1 <- hist(norm.x1, breaks = 100, plot = F)
    h2 <- hist(norm.x2, breaks = 100, plot = F) 
    a = histc(norm.x1, h1$breaks)
    b = histc(norm.x2, h1$breaks)
    c = cbind(a$cnt, b$cnt)
    d = pmin(c[,1], c[,2])
    overlap <- sum(d)
    histogram_match <- overlap/sum(a$cnt)
    print(cc); print(alpha); print(histogram_match)

    # calculate the final SPAEF
    spaef = 1- sqrt( (cc - 1)^2 + (alpha - 1)^2 + (histogram_match - 1)^2 )
    return(spaef)
}
