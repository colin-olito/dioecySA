get.maps.pathway <- function(dat, XY.YX)
{
    i <- which(XY.YX %in% colnames(dat))
    if (length(i) == 2)
    {
        ans <- dat[, c("genus", XY.YX, "root")]
    } else {
        ans <- NA
    }
    return(ans)
}

get.prop <- function(dat, XY.YX)
{
    # extract the pathway of interest
    dat.leg <- mclapply(dat, get.maps.pathway, XY.YX)
    dat.leg <- dat.leg[!is.na(dat.leg)]
    dat <- do.call(rbind, dat.leg)

    # computations in maps
    dat$sum <- dat[[XY.YX[1]]] + dat[[XY.YX[2]]]
    dat <- dat[dat$sum > 0,]
    dat$prop <- dat[[XY.YX[1]]] / dat$sum 

    return(dat)
}

get.maps <- function(genus, X)
{
    path <- paste("../analysis/analyze", X, "/", genus, sep="")

    infiles <- list.files(path, pattern="^maps_tree[0-9]{3}.csv$", 
                          full.names=T)

    # check if this genus is involved
    nruns <- length(infiles)
    if (nruns == 0)
        return(NULL)

    # combine samples from all trees
    for (i in seq_along(infiles))
    {
        dat1 <- read.csv(infiles[i], as.is=T)
        dat1$i <- NULL
        junk <- ifelse(i == 1, dat <- dat1, dat <- rbind(dat, dat1))
    }
    rm(dat1)
    dat <- cbind(genus=genus, dat)

    if (nrow(dat) != 100000)
        message(sprintf("Check if results are complete for %s (%d maps found)",
                        genus, nrow(dat)))

    return(dat)
}
