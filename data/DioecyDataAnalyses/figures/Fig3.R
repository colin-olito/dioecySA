######################################################################
### Little histograms depicting the asymmetry in transition rates. ###
######################################################################
### Asymmetry/Directionality = rXY = (qXY - qYX) / (qXY + qYX)

#--------------------------------------------------
# 1. Read in and bin results
#--------------------------------------------------

library(parallel)
options(stringsAsFactors = FALSE)

source("common-plot.R")
source("common-rates.R")
genera <- read.table("orderR.txt", as.is=T, header=F)[,1]

dat.all.4 <- mclapply(genera, get.dat, "4", mc.cores=detectCores())
names(dat.all.4) <- genera

dat.all.H <- mclapply(genera, get.dat, "H", mc.cores=detectCores())
names(dat.all.H) <- genera

dat.all.D <- mclapply(genera, get.dat, "D", mc.cores=detectCores())
names(dat.all.D) <- genera

# Pull out info for each leg alone
get.rXY <- function(dat, tnP)
{
    tnQ <- paste("q", tnP, sep="")
    rXY <- paste("r", tnP[1], sep="")

    # Keep only the columns relevant to this leg.
    # If the leg is absent, return NA instead.
    ans <- mclapply(dat, 
        function(x)
        {
            if (all(tnQ %in% names(x)))
            {
                ans <- (x[[tnQ[1]]] - x[[tnQ[2]]]) / (x[[tnQ[1]]] + x[[tnQ[2]]])
            } else {
                ans <- NA
            }
            return(ans)
        },
        mc.cores=detectCores())

    return(ans)
}

get.rate.hist <- function(dat)
{
    hist(dat, plot=F, breaks=seq(-1, 1, length.out=50))
}

rXY.hist <- list()

# 4-state

XY.YX.4 <- list(c("HG", "GH"), c("GD", "DG"), c("HM", "MH"), 
                  c("MD", "DM"), c("HD", "DH"))

for (XY.YX in XY.YX.4)
{
    rXY <- get.rXY(dat.all.4, XY.YX)
    rXY <- na.omit(unlist(rXY))
    rXY.hist[[XY.YX[1]]] <- get.rate.hist(rXY)
}

# 2-state

XY.YX.2 <- list(c("HO", "OH"), c("OD", "DO"))

XY.YX <- XY.YX.2[[1]]
rXY <- get.rXY(dat.all.H, XY.YX)
rXY <- na.omit(unlist(rXY))
rXY.hist[[XY.YX[1]]] <- get.rate.hist(rXY)

XY.YX <- XY.YX.2[[2]]
rXY <- get.rXY(dat.all.D, XY.YX)
rXY <- na.omit(unlist(rXY))
rXY.hist[[XY.YX[1]]] <- get.rate.hist(rXY)

#--------------------------------------------------
# 2. Write out results for use in Fig5
#--------------------------------------------------

library(coda)
get.ci.width <- function(x)
{
    if (is.na(x[1]))
    {
        ans <- NA
    } else {
        ans <- HPDinterval(as.mcmc(x), prob=0.90)
        ans <- ans[2] - ans[1]
    }
    return(ans)
}

XY.YX.all <- c(XY.YX.4, XY.YX.2)

rpXY <- data.frame(matrix(NA, nrow=length(genera), ncol=length(XY.YX.all)*2))
rownames(rpXY) <- genera
XY <- sapply(XY.YX.all, function(x) x[1])
colnames(rpXY) <- c(XY, paste(XY, ".90", sep=""))

# Median, and width of 90% CI

# 4-state
for (XY.YX in XY.YX.4)
{
    rXY <- get.rXY(dat.all.4, XY.YX)
    ans.med <- sapply(rXY, median)
    ans.ci <- sapply(rXY, get.ci.width)
    rpXY[names(ans.med), XY.YX[1]] <- ans.med
    rpXY[names(ans.ci), paste(XY.YX[1], ".90", sep="")] <- ans.ci
}

# HO
XY.YX <- XY.YX.2[[1]]
rXY <- get.rXY(dat.all.H, XY.YX)
ans.med <- sapply(rXY, median)
ans.ci <- sapply(rXY, get.ci.width)
rpXY[names(ans.med), XY.YX[1]] <- ans.med
rpXY[names(ans.ci), paste(XY.YX[1], ".90", sep="")] <- ans.ci

# OD
XY.YX <- XY.YX.2[[2]]
rXY <- get.rXY(dat.all.D, XY.YX)
ans.med <- sapply(rXY, median)
ans.ci <- sapply(rXY, get.ci.width)
rpXY[names(ans.med), XY.YX[1]] <- ans.med
rpXY[names(ans.ci), paste(XY.YX[1], ".90", sep="")] <- ans.ci

write.table(rpXY, file="Fig5a.dat")

#--------------------------------------------------
# 3. Make the little histograms
#--------------------------------------------------

XY.YX.all <- c(XY.YX.4, XY.YX.2)

ylim <- c(0, max(sapply(rXY.hist, function(x) max(x$density))))

pdf("Fig3.pdf", width=6, height=5)
par(mar = c(0,0,1,0), oma=c(0,0,0,0))
for (XY.YX in XY.YX.all)
{
    XY <- XY.YX[1]
    ans <- rXY.hist[[XY]]
    plot(ans, freq = F, col="black", axes=F, xlab="", ylab="", 
         main=XY, ylim=ylim)
}
dev.off()
