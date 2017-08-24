### For each pathway,
###     show each rate estimate
###     show differences in rate estimates; highlight when significant

library(parallel)
source("common-plot.R")
source("common-rates.R")

################
### 4-states ###
################

genera <- read.table("orderR.txt", as.is=T, header=F)[,1]

#--------------------------------------------------
# Assemble data
#--------------------------------------------------

# forwards, backwards, difference
dat.all <- mclapply(genera, get.dat, "4", mc.cores=detectCores())
names(dat.all) <- genera

# prior
prior.all <- mclapply(genera, get.prior, "4", mc.cores=detectCores())
names(prior.all) <- genera

# That was for all the pathways and legs together.
# Now pull out info for each leg alone.
get.leg <- function(dat, tnP)
{
    tnQ <- paste("q", tnP, sep="")
    XP <- paste(tnP, collapse=".")
    XQ <- paste("q", XP, sep="")

    # Keep only the columns relevant to this leg.
    # If the leg is absent, return NA instead.
    datL <- mclapply(dat, 
        function(x)
        {
            if (all(tnQ %in% names(x)))
            {
                ans <- x[,tnQ]
                ans[[XQ]] <- x[[tnQ[1]]] - x[[tnQ[2]]]
            } else {
                ans <- NA
            }
            return(ans)
        },
        mc.cores=detectCores())

    # Test if rates are significantly greater than zero, or different ###
    sigL <- data.frame(unlist(mclapply(datL, function(x) get.ci(x, tnQ[1], p)[1] > eps, mc.cores=detectCores())), 
                       unlist(mclapply(datL, function(x) get.ci(x, tnQ[2], p)[1] > eps, mc.cores=detectCores()))
                          )
    names(sigL) <- tnP
    rownames(sigL) <- genP
    sigL[[XP]] <- unlist(mclapply(datL, 
        function(x)
        {
            if (!is.na(x))
                get.sig.diff(data.frame(d = x[[tnQ[1]]] - x[[tnQ[2]]]), "d", p)
            else
                NA
        },
        mc.cores=detectCores()))

    return(list(datL=datL, sigL=sigL))
}

#--------------------------------------------------
# HD pathway
#--------------------------------------------------

### Extract the data

# Get genera involved in this pathway

path.names <- c("qHD", "qDH")
genP <- read.table("orderHD.txt", as.is=T, header=F)[,1]
genP <- genP[which(genP %in% genera)]

# Note which genera are large

font.names <- rep(3, length(genP))
font.names[which(genP %in% genera.large)] <- 4

# Filter the data down to the relevant genera

i <- sapply(genP, function(x) which(names(dat.all) == x))
datP <- dat.all[i]
priorP <- prior.all[genP]

# Get the data for each leg

leg1 <- get.leg(datP, c("HD", "DH"))
leg1$sigL[,1:2] <- 0  # use significance shading only for rate difference

### Plot the figure

plotfile <- "FigS9.pdf"
pdf(file=plotfile, width=4, height=length(genP)/2.75, pointsize=16)
par(mfcol=c(length(genP), 3), mar=c(0.6,0,0.1,1), oma=c(1,7,3,0))

plot.set(leg1$datL, priorP, leg1$sigL, "HD", with.names=T, font.names=font.names)
mtext(expression(italic(q[HD])), side=3, line=0, cex=0.8, outer=T, adj=0.1)
plot.set(leg1$datL, priorP, leg1$sigL, "DH")
mtext(expression(italic(q[DH])), side=3, line=0, cex=0.8, outer=T, adj=0.45)
plot.set(leg1$datL, priorP, leg1$sigL, "HD.DH")
mtext(expression(italic(q[HD] - q[DH])), side=3, line=0, cex=0.8, outer=T, adj=0.9)

mtext("Direct pathway", outer=T, line=2, cex=0.75, font=2, adj=0.5)
dev.off()

#--------------------------------------------------
# HMD pathway
#--------------------------------------------------

### Extract the data

# Get genera involved in either leg of this pathway
path.names <- c("qHM", "qMH", "qMD", "qDM")
genP <- read.table("orderHMD.txt", as.is=T, header=F)[,1]
genP <- genP[which(genP %in% genera)]

# Note which genera are large

font.names <- rep(3, length(genP))
font.names[which(genP %in% genera.large)] <- 4

# Filter the data down to the relevant genera

i <- sapply(genP, function(x) which(names(dat.all) == x))
datP <- dat.all[i]
priorP <- prior.all[genP]

# Get the data for each leg

leg1 <- get.leg(datP, c("HM", "MH"))
leg2 <- get.leg(datP, c("MD", "DM"))
# use significance shading only for rate difference
leg1$sigL[,1:2] <- 0
leg2$sigL[,1:2] <- 0

### Plot the figure

plotfile <- "FigS8.pdf"
pdf(file=plotfile, width=6, height=length(genP)/2.75, pointsize=16)
par(mfcol=c(length(genP), 6), mar=c(0.6,0,0.1,1), oma=c(1,7,3.5,0))

plot.set(leg1$datL, priorP, leg1$sigL, "HM", with.names=T, font.names=font.names)
mtext(expression(italic(q[HM])), side=3, line=0.4, cex=0.8, outer=T, adj=0.01)
plot.set(leg1$datL, priorP, leg1$sigL, "MH")
mtext(expression(italic(q[MH])), side=3, line=0.4, cex=0.8, outer=T, adj=0.2)
plot.set(leg1$datL, priorP, leg1$sigL, "HM.MH")
mtext(expression(italic(q[HM] - q[MH])), side=3, line=0.4, cex=0.8, outer=T, adj=0.38)

plot.set(leg2$datL, priorP, leg2$sigL, "MD")
mtext(expression(italic(q[MD])), side=3, line=0.4, cex=0.8, outer=T, adj=0.57)
plot.set(leg2$datL, priorP, leg2$sigL, "DM")
mtext(expression(italic(q[DM])), side=3, line=0.4, cex=0.8, outer=T, adj=0.74)
plot.set(leg2$datL, priorP, leg2$sigL, "MD.DM")
mtext(expression(italic(q[MD] - q[DM])), side=3, line=0.4, cex=0.8, outer=T, adj=0.98)

mtext("Monomorphic pathway", outer=T, line=2.5, cex=0.75, font=2, adj=0.5)
dev.off()

#--------------------------------------------------
# HGD pathway
#--------------------------------------------------

### Extract the data

# Get genera involved in either leg of this pathway

path.names <- c("qHG", "qGH", "qGD", "qDG")
genP <- read.table("orderHGD.txt", as.is=T, header=F)[,1]
genP <- genP[which(genP %in% genera)]

# Note which genera are large

font.names <- rep(3, length(genP))
font.names[which(genP %in% genera.large)] <- 4

# Filter the data down to the relevant genera

i <- sapply(genP, function(x) which(names(dat.all) == x))
datP <- dat.all[i]
priorP <- prior.all[genP]

# Get the data for each leg

leg1 <- get.leg(datP, c("HG", "GH"))
leg2 <- get.leg(datP, c("GD", "DG"))
# use significance shading only for rate difference
leg1$sigL[,1:2] <- 0
leg2$sigL[,1:2] <- 0

### Plot the figure

plotfile <- "FigS7.pdf"
pdf(file=plotfile, width=6, height=length(genP)/2.75, pointsize=16)
par(mfcol=c(length(genP), 6), mar=c(0.6,0,0.1,1), oma=c(1,7,3.5,0))

plot.set(leg1$datL, priorP, leg1$sigL, "HG", with.names=T, font.names=font.names)
mtext(expression(italic(q[HG])), side=3, line=0.4, cex=0.8, outer=T, adj=0.01)
plot.set(leg1$datL, priorP, leg1$sigL, "GH")
mtext(expression(italic(q[GH])), side=3, line=0.4, cex=0.8, outer=T, adj=0.2)
plot.set(leg1$datL, priorP, leg1$sigL, "HG.GH")
mtext(expression(italic(q[HG] - q[GH])), side=3, line=0.4, cex=0.8, outer=T, adj=0.38)

plot.set(leg2$datL, priorP, leg2$sigL, "GD")
mtext(expression(italic(q[GD])), side=3, line=0.4, cex=0.8, outer=T, adj=0.57)
plot.set(leg2$datL, priorP, leg2$sigL, "DG")
mtext(expression(italic(q[DG])), side=3, line=0.4, cex=0.8, outer=T, adj=0.74)
plot.set(leg2$datL, priorP, leg2$sigL, "GD.DG")
mtext(expression(italic(q[GD] - q[DG])), side=3, line=0.4, cex=0.8, outer=T, adj=0.98)

mtext("Dimorphic pathway", outer=T, line=2.5, cex=0.75, font=2, adj=0.5)
dev.off()

##########
### HO ###
##########

genera <- read.table("orderH.txt", as.is=T, header=F)[,1]

# Note which genera are large
font.names <- rep(3, length(genera))
font.names[which(genera %in% genera.large)] <- 4

#--------------------------------------------------
# Assemble data
#--------------------------------------------------

# forwards, backwards, difference
dat.all <- mclapply(genera, get.dat, "H", mc.cores=detectCores())
names(dat.all) <- genera

# prior
prior.all <- mclapply(genera, get.prior, "H", mc.cores=detectCores())
names(prior.all) <- genera

### Test if rates are significantly greater than zero, or different ###

sig.all <- data.frame("HO" = unlist(mclapply(dat.all, 
                             function(x) get.ci(x, "qHO", p)[1] > eps,
                             mc.cores=detectCores())), 
                      "OH" = unlist(mclapply(dat.all, 
                             function(x) get.ci(x, "qOH", p)[1] > eps,
                             mc.cores = detectCores()))
                      )
sig.all[["HO.OH"]] <- unlist(mclapply(dat.all, get.sig.diff, "qHO.OH", p, mc.cores=detectCores()))
rownames(sig.all) <- genera

# Use significance shading only for rate difference
sig.all[,1:2] <- FALSE

#--------------------------------------------------
# Draw the figure
#--------------------------------------------------

plotfile <- "FigS5.pdf"
pdf(file=plotfile, width=4, height=length(genera)/2.75, pointsize=16)
par(mfcol=c(length(genera), 3), mar=c(0.6,0,0.1,1), oma=c(1,6.2,3.5,0))

plot.set(dat.all, prior.all, sig.all, "HO", with.names=T, font.names=font.names)
mtext(expression(italic(q[HO])), side=3, line=0.4, cex=0.8, outer=T, adj=0.1)

plot.set(dat.all, prior.all, sig.all, "OH")
mtext(expression(italic(q[OH])), side=3, line=0.4, cex=0.8, outer=T, adj=0.45)

plot.set(dat.all, prior.all, sig.all, "HO.OH")
mtext(expression(italic(q[HO] - q[OH])), side=3, line=0.4, cex=0.8, outer=T, adj=0.9)

mtext("Transition rates from/to hermaphroditism", outer=T, line=2, cex=0.75, font=2, adj=1.5)
dev.off()

##########
### OD ###
##########

genera <- read.table("orderD.txt", as.is=T, header=F)[,1]

# Note which genera are large
font.names <- rep(3, length(genera))
font.names[which(genera %in% genera.large)] <- 4

#--------------------------------------------------
# Assemble data
#--------------------------------------------------

# forwards, backwards, difference
dat.all <- mclapply(genera, get.dat, "D", mc.cores=detectCores())
names(dat.all) <- genera

# prior
prior.all <- mclapply(genera, get.prior, "D", mc.cores=detectCores())
names(prior.all) <- genera

### Test if rates are significantly greater than zero, or different ###

sig.all <- data.frame("OD" = unlist(mclapply(dat.all, 
                             function(x) get.ci(x, "qOD", p)[1] > eps,
                             mc.cores=detectCores())), 
                      "DO" = unlist(mclapply(dat.all, 
                             function(x) get.ci(x, "qDO", p)[1] > eps,
                             mc.cores=detectCores()))
                      )
sig.all[["OD.DO"]] <- unlist(mclapply(dat.all, get.sig.diff, "qOD.DO", p, mc.cores=detectCores()))
rownames(sig.all) <- genera

# Use significance shading only for rate difference
sig.all[,1:2] <- FALSE

#--------------------------------------------------
# Draw the figure
#--------------------------------------------------

plotfile <- "FigS6.pdf"
pdf(file=plotfile, width=4, height=length(genera)/2.75, pointsize=16)
par(mfcol=c(length(genera), 3), mar=c(0.6,0,0.1,1), oma=c(1.5,7,3.5,0))

plot.set(dat.all, prior.all, sig.all, "OD", with.names=T, font.names=font.names)
mtext(expression(italic(q[OD])), side=3, line=0.4, cex=0.8, outer=T, adj=0.1)

plot.set(dat.all, prior.all, sig.all, "DO")
mtext(expression(italic(q[DO])), side=3, line=0.4, cex=0.8, outer=T, adj=0.45)

plot.set(dat.all, prior.all, sig.all, "OD.DO")
mtext(expression(italic(q[OD] - q[DO])), side=3, line=0.4, cex=0.8, outer=T, adj=0.9)

mtext("Transition rates to/from dioecy", outer=T, line=2, cex=0.75, font=2, adj=0)
dev.off()
