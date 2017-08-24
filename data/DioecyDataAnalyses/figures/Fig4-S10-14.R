#####################################################################
### Little histograms depicting the asymmetry in each transition. ###
#####################################################################
### Asymmetry/Directionality = dXY = nXY / (nXY + nYX)
### Main results: all genera combined
### Supp results: each genus separately (without symmetric rates)

#--------------------------------------------------
# 1. Read in and bin results
#    For all transitions
#--------------------------------------------------

library(parallel)
source("common-plot.R")
source("common-maps.R")
options(stringsAsFactors = FALSE)
genera <- read.table("orderR.txt", as.is=T, header=F)[,1]

get.prop.hist <- function(dat)
{
    hist(dat$prop, plot=F, breaks=seq(0, 1, length.out=30))
}

### 4-state mappings ###

XY.YX.4 <- list(c("HG", "GH"), c("GD", "DG"), c("HM", "MH"), 
                  c("MD", "DM"), c("HD", "DH"))

dat.real <- mclapply(genera, get.maps, "4", mc.cores=detectCores())
dat.prop.real <- list()
for (XY.YX in XY.YX.4)
{
    dat.prop.real[[XY.YX[1]]] <- get.prop(dat.real, XY.YX)
}

rm(dat.real)

### 2-state mappings ###

XY.YX.2 <- list(c("HO", "OH"), c("OD", "DO"))


dat.H <- mclapply(genera, get.maps, "H", mc.cores=detectCores())
dat.D <- mclapply(genera, get.maps, "D", mc.cores=detectCores())

dat.prop.real.H <- get.prop(dat.H, XY.YX.2[[1]])
dat.prop.real.D <- get.prop(dat.D, XY.YX.2[[2]])

### 4 and 2 together ###
dat.prop.real <- c(dat.prop.real, 
                   list("HO"=dat.prop.real.H), 
                   list("OD"=dat.prop.real.D))
hist.real <- lapply(dat.prop.real, get.prop.hist)
XY.YX.all <- c(XY.YX.4, XY.YX.2)

rm(dat.H, dat.D, dat.prop.real.H, dat.prop.real.D)

#--------------------------------------------------
# 2. Plot dXY results for all genera combined
#--------------------------------------------------

ymax <- max(sapply(hist.real, function(x) max(x$density)))

pdf(file="Fig4.pdf", width=3, height=2)
for (XY.YX in XY.YX.all)
{
    X <- substr(XY.YX[1], 1, 1)
    Y <- substr(XY.YX[1], 2, 2)
    ans <- hist.real[[XY.YX[1]]]

    par(mar = c(0,0,0,0), oma=c(0,0,0,0))
    plot(ans, freq = F, col="black", axes=F, xlab="", ylab="", main="", 
         ylim=c(0, ymax*1.12))
    text(1, ymax, paste(X, Y, sep = " -> "), adj=c(1, -1))
    text(0, ymax, paste(Y, X, sep = " -> "), adj=c(0, -1))
}
dev.off()

#--------------------------------------------------
# 3. Write out results for use in Fig5b
#--------------------------------------------------
# pXY = mean(dXY)

library(coda)
get.ci.width <- function(x, g)
{
    x <- subset(x, genus==g)$prop
    if (is.na(x[1]))
    {
        ans <- NA
    } else {
        ans <- HPDinterval(as.mcmc(x), prob=0.90)
        ans <- ans[2] - ans[1]
    }
    return(ans)
}

out <- data.frame(matrix(NA, nrow=length(genera), ncol=14))
row.names(out) <- genera

XY <- sapply(XY.YX.all, function(x) x[1])
colnames(out) <- c(XY, paste(XY, ".90", sep=""))

for (genus in genera)
{
    ans.med <- sapply(dat.prop.real, function(x, g) 
                          mean(subset(x, genus==g)$prop), genus)
    ans.ci <- sapply(dat.prop.real, function(x, g) get.ci.width(x, g), genus)
    out[genus,] <- c(ans.med, ans.ci)
}

write.table(out, file="Fig5b.dat")

#--------------------------------------------------
# 4. Plot each genus separately
#    Real rates only
#--------------------------------------------------

pdf(file="FigS10-14.pdf", width=12, height=8)
for (XY2 in XY.YX.all)
{
    XY <- XY2[1]
    dat1 <- dat.prop.real[[XY]]

    par(mfrow=c(6,7), mar=c(1,1,1,0), oma=c(0,0,1.2,0))

    for (gen in sort(unique(dat1$genus)))
    {
        dat.gen.real <- subset(dat1, genus == gen)
        ans.real <- get.prop.hist(dat.gen.real)

        # ylim <- c(0, max(ans.real$density))
        ylim <- c(0, 26)
        plot(ans.real, col="black", ylim=ylim, main="", freq=F, axes=F)
        title(gen, font.main=3)

        n.maps <- round(nrow(dat.gen.real) / 100000, 2)
        if (n.maps != 1)
            text(0.2, ylim[2]*0.95, n.maps)
    }
    mtext(XY, line=0, outer=T)

}
dev.off()
