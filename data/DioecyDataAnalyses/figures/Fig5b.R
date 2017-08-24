############################################################
### One point per genus per state pair                   ###
### Average proportion of state changes that are forward ###
############################################################
# pXY = mean(dXY)

#--------------------------------------------------
# 1. Read in results
#--------------------------------------------------

options(stringsAsFactors = FALSE)

dat <- read.table("Fig5b.dat") # computed in Fig4-S10-14.R
genera <- rownames(dat)

library(viridis)
col.genera <- rev(viridis(nrow(dat)))
names(col.genera) <- genera

library(plotrix)

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

XY <- rev(c("HG", "HM", "HD", "GD", "MD", "HO", "OD"))

# Separate median from CI info
dat.conf <- dat[, paste(XY, ".90", sep="")] <= 0.25
dat <- dat[, XY]
colnames(dat.conf) <- colnames(dat)

#--------------------------------------------------
# 2. Make the plot
#--------------------------------------------------

pdf(file="Fig5b.pdf", width=9, height=4.5)
par(xpd = NA, oma=c(0,0,0,3.5), mar=c(3.5,2,0.1,2))

plot(NA, xlim=c(0, 1), ylim=c(0.9, 7.1), xlab="", ylab="", axes=F)
axis(side=2, labels=colnames(dat), at=1:7)
axis(side=1)
mtext("proportion of transitions", 1, line=2.5)
box()
segments(0.5, 0.7, 0.5, ncol(dat)+0.35, lty=3)

text.font <- structure(rep(3, length(genera)), names=genera) 

for (i in seq(ncol(dat)))
{
    x1 <- dat[, i]
    set.seed(3)
    y1 <- rep(i, length(x1)) + 0.12 + rnorm(length(x1), 0, 0.05)
    y2 <- i - 0.12

    cex <- structure(rep(1, length(genera)), names=genera)
    g <- names(which(dat.conf[, XY[i]]))
    cex[g] <- 1.4
    text.font[g] <- 4

    points(x1, y1, col=col.genera, lwd=2, cex=cex)

    se <- stderr(x1)
    plotCI(mean(x1, na.rm=T), y2, se, err="x", add=TRUE, pch=16, cex=1.4)
}

legend(par("usr")[2], par("usr")[4]+0.15, genera, text.col=col.genera, 
       text.font=text.font, cex=0.56, bty="n")

dev.off()
