######################################################
# Barplots for each genus, showing state frequencies #
######################################################

options(stringsAsFactors = FALSE)
genera <- read.table("orderR.txt", as.is=T, header=F)[,1]

source("common-plot.R")
cols <- col.state

#--------------------------------------------------
# Read in the results
#--------------------------------------------------

genera <- rev(genera)
tnames <- c("H", "M", "G", "D")

### Root

infile <- "../analysis/analyze4/summaries/root_states.csv"
dat <- read.csv(infile, as.is=T)
rownames(dat) <- dat$X
dat$X <- NULL
dat.r <- dat[genera, tnames]
rm(dat)
dat.r <- t(as.matrix(dat.r))
dat.r[which(is.na(dat.r))] <- 0

sig.r <- (dat.r > p)

### Tips

infile <- "../analysis/analyze4/summaries/tip_states.csv"
dat <- read.csv(infile, as.is=T)
rownames(dat) <- dat$X
dat$X <- NULL
dat.t <- dat[genera, tnames]
rm(dat)
dat.t <- t(as.matrix(dat.t))
dat.t <- apply(dat.t, 2, function(x) x / sum(x))

### Durations

infile <- "../analysis/analyze4/summaries/state_durations.csv"
dat <- read.csv(infile, as.is=T)
rownames(dat) <- dat$X
dat$X <- NULL
dat.d <- dat[genera, tnames]
rm(dat)
dat.d <- t(as.matrix(dat.d))
dat.d[which(is.na(dat.d))] <- 0

### Numbers of tips
infile <- "ntips.csv"
dat <- read.csv(infile, as.is=T)
dat.n <- structure(dat$tips.num, names=dat$Genus)

#--------------------------------------------------
# Make the plots
#--------------------------------------------------

draw.rect <- function(X)
{
    y <- c(NA, NA)
    i <- which(sig.r[X,])
    if (length(i) > 0)
    {
        y[1] <- min(i) * 2 - 1.5
        y[2] <- max(i) * 2 + 0.5
        rect(xlim[1], y[1], xlim[2], y[2], col=col.root[X],
             border=col.root[X])
    }
}

pdf(file="FigS1.pdf", width=5.4, height=length(genera)/2.1, pointsize=16)
par(mar=c(0,0,0,0), oma=c(0,0,0,1))

xlim <- c(-0.1, 1.1)

vw <- c(0, 0.31, 0.1, 0.19, 0.19, 0.19)
screen.num <- split.screen(matrix(c(cumsum(vw)[-6], cumsum(vw)[-1], rep(0,5), rep(1,5)), ncol=4))

# highlight large genera
font.gen <- rep(3, length(genera))
font.gen[which(genera %in% genera.large)] <- 4
vr <- seq(0.005, 0.992, length.out=length(genera))
screen(1)
text(1, vr, genera, pos=2, font=font.gen, cex=1)

# number of species per genus
screen(2)
text(1, vr, dat.n[genera], pos=2)

# tip states
screen(3)
barplot(dat.t, axes=F, col=cols[rownames(dat.t)], horiz=T, space=1, las=3,
        font.axis=3, cex.names=1.3, xlim=xlim, 
        names.arg=rep("", ncol(dat.t)))
title("Tips", line=-1.5, font.main=2, cex.main=1)

# root states
screen(4)
barplot(dat.r, axes=F, col=cols[rownames(dat.r)], horiz=T, space=1, las=3,
        font.axis=3, cex.names=1.3, xlim=xlim, 
        names.arg=rep("", ncol(dat.t)))
draw.rect("H")
draw.rect("M")
draw.rect("G")
draw.rect("D")
barplot(dat.r, axes=F, col=cols[rownames(dat.r)], horiz=T, space=1, las=2,
        font.axis=3, cex.names=1.3, add=T, names.arg=rep("", ncol(dat.t)))
title("Root", line=-1.5, font.main=2, cex.main=1)

# durations in each state
screen(5)
barplot(dat.d, axes=F, col=cols[rownames(dat.d)], horiz=T, space=1, las=3,
        font.axis=3, cex.names=1.3, xlim=xlim, 
        names.arg=rep("", ncol(dat.d)))
title("Duration", line=-1.5, font.main=2, cex.main=1)

close.screen(all = TRUE)
dev.off()
