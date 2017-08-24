######################################################
### Barplot for number of transitions of each type ###
######################################################
# nXY over all mappings, scaled by total number of transitions

library(parallel)
library(plotrix)

source("common-plot.R", chdir=T)
source("common-maps.R")
options(stringsAsFactors = FALSE)
genera <- read.table("orderR.txt", as.is=T, header=F)[,1]

XY.YX.4 <- list(c("HG", "GH"), c("GD", "DG"), c("HM", "MH"), 
                  c("MD", "DM"), c("HD", "DH"), c("MG", "GM"))
XY.YX.2 <- list(c("HO", "OH"), c("OD", "DO"))
XY.4 <- sapply(XY.YX.4, function(x) x[1])   # forward
XY.2 <- sapply(XY.YX.2, function(x) x[1])   # forward

#--------------------------------------------------
# 1. Prepare the summary data file
#--------------------------------------------------

dat.4 <- mclapply(genera, get.maps, "4", mc.cores=detectCores())
dat.H <- mclapply(genera, get.maps, "H", mc.cores=detectCores())
dat.D <- mclapply(genera, get.maps, "D", mc.cores=detectCores())

countXY <- function(XY, dat)
{
    func <- function(dat.gen, XY)
        ifelse (XY %in% names(dat.gen), sum(dat.gen[,XY]), 0)
    sum(sapply(dat, func, XY))
}

### All genera together

ans4 <- sapply(unlist(XY.YX.4), countXY, dat.4)
ansH <- sapply(XY.YX.2[[1]], countXY, dat.H)
ansD <- sapply(XY.YX.2[[2]], countXY, dat.D)

write.table(c(ans4, ansH, ansD), file="total_changes.csv", col.names=F, sep=",")

#--------------------------------------------------
# 2. Plot the results
#--------------------------------------------------

dat0 <- read.csv("total_changes.csv", header=F)
dat0 <- structure(dat0[,2], names=dat0[,1])

nmaps <- 100000 * length(genera)
dat1 <- dat0 / nmaps
dat <- matrix(dat1, nrow=2, 
              dimnames = list(c("more SD", "less SD"), c(XY.4, XY.2)))

pdf(file="Fig2.pdf", height=3, width=4.5)
par(mar=c(2,4,1,0.5))
barplot(dat, legend.text=T, args.legend=list(x="topleft"), 
        ylab="average number of transitions per mapping", 
        col=c("#888888", "#eeeeee"))
dev.off()
