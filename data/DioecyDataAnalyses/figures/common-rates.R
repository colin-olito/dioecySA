col.prior <- mygray # or "black"

# Choose the threshold for "greater than 0"
eps <- 1e-3

#--------------------------------------------------
# Assemble data
#--------------------------------------------------

get.dat <- function(genus, X)
{
    # check if results exist
    path <- paste("../analysis/analyze", X, "/", genus, sep="")
    infiles <- list.files(path, pattern="^mcmc_tree[0-9]{3}.csv$", full.names=T)
    nruns <- length(infiles)
    if (nruns == 0)
    {
        message(paste(genus, ": ", X, " rates files not found", sep=""))
        return(NULL)
    }

    # combine samples from all trees
    for (i in seq_along(infiles))
    {
        dat1 <- read.csv(infiles[i])
        dat1$i <- dat1$p <- NULL
        junk <- ifelse(i == 1, dat <- dat1, dat <- rbind(dat, dat1))
    }

    # compute rate difference
    if (X == "H")
    {
        dat$qHO.OH <- dat$qHO - dat$qOH
    } else if (X == "D") {
        dat$qOD.DO <- dat$qOD - dat$qDO
    }

    # discard speciation and extinction rates
    i <- which(names(dat) %in% c("lambda", "mu"))
    if (length(i) > 0)
        dat <- dat[,-i]

    if (nrow(dat) != 100000)
        message(paste("Check if results are complete for", genus))

    return(dat)
}

get.prior <- function(genus, X)
{
    priors <- c("4" = 0.35, "H" = 0.6, "D" = 0.4)
    prior.mean <- 1/priors[X]
    return(prior.mean)
}

# could use coda instead
get.ci <- function(dat, varname, p)
{
    if (!is.na(dat))
    {
        cis <- diversitree:::hdr(dat[[varname]], p=p)
        names(cis) <- paste(p, c("Low", "Hi"), sep="")
    } else {
        cis <- NA
    }
    return(cis)
}

get.sig.diff <- function(dat, varname, p)
{
    dat <- dat[[varname]]
    p0 <- sum(dat > 0) / length(dat)

    if (p0 > p)
    {
        ans <- 1
    } else if (p0 < 1-p) {
        ans <- 2
    } else {
        ans <- 0
    }

    return(ans)    
}

#--------------------------------------------------
# Plotting functions
#--------------------------------------------------

### Plot one histogram of one rate ###

# (If sig is not NA, it should be the color of the background box.)
plot.hist <- function(dat, prior.mean, col="black", sig=NA)
{
    # Using hist rather than dens to avoid downturn at 0.
    dat.hist <- hist(dat, breaks=100, plot=F)
    x <- dat.hist$mids
    y <- dat.hist$density
    xlim <- c(0, 3)     # some multiple of root age
    ylim <- c(0, max(y))
    i <- which(x <= xlim[2])
    x <- x[i]
    y <- y[i]

    plot(x, y, xlim=xlim, ylim=ylim, axes=F, ann=F, type="l")
    axis(side=1, at=seq(0, xlim[2], by=0.5), label=F, tcl=-0.19)
    axis(side=1, at=seq(0, xlim[2]), cex.axis=0.6, tcl=-0.19, padj=-3.5)
    if (!is.na(sig))
    {
        rect(xlim[1], ylim[1], xlim[2], ylim[2], col=sig, border=NA)
    }
    polygon(c(0, 0, x, max(x)), 
            c(0, y[1], y, 0), 
            col=col, border=col)

    x.p <- seq(xlim[1], xlim[2], length.out=1000)
    lines(x.p, dexp(x.p, 1/prior.mean), col=col.prior)
}

### Plot one histogram of rate difference ###

plot.hist.diff <- function(dat, prior.mean, sig=NA)
{
    # density looks smoother than hist at 0
    dat.dens <- density(dat)
    x <- dat.dens$x
    y <- dat.dens$y
    xlim <- diversitree:::hdr(dat, p=p)
    ylim <- c(0, max(y))
    i <- which(x >= xlim[1] & x <= xlim[2])
    x <- x[i]
    y <- y[i]

    plot(x, y, xlim=xlim, ylim=ylim, axes=F, ann=F, type="l")

    at.all <- seq(-10, 10, by=0.5)
    at.int <- seq(-10, 10, by=1)
    label.int <- as.character(at.int)
    i <- which(at.int >= xlim[1] & at.int <= xlim[2])
    axis(side=1, at=at.all, label=F, tcl=-0.19)
    mtext(side=1, at=at.int[i], text=label.int[i], cex=0.4, line=-0.2)

    if (!is.na(sig))
    {
        rect(xlim[1], ylim[1], xlim[2], ylim[2], col=sig, border=NA)
    }

    i <- which(x >= 0)
    if (length(i) > 0)
    {
        x0 <- max(min(x[i]), 0)  # in case hist doesn't include 0
        polygon(c(x0, x0, x[i], max(x[i])), 
                c(0, y[i[1]], y[i], 0), 
                col=col.dark[1], border=col.dark[1])
    }
    i <- which(x <= 0)
    if (length(i) > 0)
    {
        x0 <- min(max(x[i]), 0)  # in case hist doesn't include 0
        polygon(c(x0, min(x[i]), x[i], x0), 
                c(0, 0, y[i], y[i[length(i)]]), 
                col=col.dark[2], border=col.dark[2])
    }

    x.p <- seq(0, max(abs(xlim)), length.out=1000)
    i <- which(x.p < xlim[2])
    lines(x.p[i], dexp(x.p[i], 1/prior.mean), col=col.prior)
    i <- which(-x.p > xlim[1])
    lines(-x.p[i], dexp(x.p[i], 1/prior.mean), col=col.prior)
}

### Assemble all histograms ###

# dat = list of dataframes
# prior = list of scalars
plot.set <- function(dat, prior, sig.all, XY, screen.num=NA, 
                     with.names=FALSE, font.names=3)
{
    if (XY %in% c("HO", "OD", "HD", "HM", "MD", "HG", "GD"))
    {
        # for the forwards columns
        col <- col.dark[1]
        sig.col <- col.light[1]
    } else if (XY %in% c("OH", "DO", "DH", "MH", "DM", "GH", "DG")) {
        # for the backwards columns
        col <- col.dark[2]
        sig.col <- col.light[2]
    } else {
        # for the rate difference column
        col <- NA
        sig.col <- col.light
    }
    if (with.names & length(font.names) == 1)
    {
        font.names <- rep(font.names, length(dat))
    }

    for (i in seq_along(dat))
    {
        genus <- names(dat)[i]
        dat1 <- dat[[i]]
        if (any(!is.na(dat1)))
        {
            prior.mean <- prior[[i]]
            sig <- sig.all[genus, XY]
            if (sig & XY %in% c("HO.OH", "OD.DO", "HD.DH", "HM.MH", "MD.DM",
                                "HG.GH", "GD.DG"))
            {
                # rate difference
                sig <- sig.col[sig]
            } else if (sig) {
                # single rate
                sig <- sig.col
            }

            q.names <- colnames(dat1)[grep("q", colnames(dat1))]
            j <- which(q.names == paste("q", XY, sep=""))

            if (!is.na(screen.num))
            {
                screen(screen.num)
                par(mar=c(0.1,0,0.1,1))
                screen.num <- screen.num + 1
            }

            if (XY %in% c("HO.OH", "OD.DO", "HD.DH", "HM.MH", "MD.DM",
                          "HG.GH", "HG.GD", "GD.DG"))
            {
                plot.hist.diff(dat1[,j], prior.mean, sig)
            } else {
                plot.hist(dat1[,j], prior.mean, col, sig)
            }
        } else {
            if (is.na(screen.num))
                frame()
        }

        if (with.names)
            mtext(genus, 2, font=font.names[i], line=1.5, las=2, cex=0.7)
    }
}
