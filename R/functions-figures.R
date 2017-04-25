###############
# DEPENDENCIES
###############
library(plyr)
library(plotrix)

#######################
# AUXILLIARY FUNCTIONS
#######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('output/figures', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}



####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
plotGrid  <-  function(lineCol='white',...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
}


#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}


#' Creates a vector of sf values spaced across a slice of funnelplots for a given sm value
#'
#' @title Creates a vector of sf values spaced across a slice of funnelplots for a given sm value
#' @param sm, selection coefficient via male function.
#' @param C, population selfing rate.
#' @param delta, population inbreeding depression.
#' @param h, dominance coefficients. Default is additive fitness effects (hf = hm = 1/2). 
#'           Can also provide arbitrary value h for the case of dominance reversal (h = hf = hm < 1/2)
#' @return A numeric vector.
#' @author Colin Olito.
#' @export
funnelSlice  <-  function(sm, C, delta, h = 0.5) { 
    if(h == 0.5){
        if(Inv.a.add(sm, C, delta) >= 1) {
        y  <-  c(Inv.A.add(sm, C, delta) - 0.1*(1 - Inv.A.add(sm, C, delta)),
                     Inv.A.add(sm, C, delta) +   (1 - Inv.A.add(sm, C, delta))/10,
                     Inv.A.add(sm, C, delta) +   (1 - Inv.A.add(sm, C, delta))/5,
                     Inv.A.add(sm, C, delta) +   (1 - Inv.A.add(sm, C, delta))/2,
                     1
                    )
            }
        if(Inv.a.add(sm, C, delta) < 1) {
            y  <-  c(Inv.A.add(sm, C, delta) - 0.1*(Inv.a.add(sm, C, delta) - Inv.A.add(sm, C, delta)),
                     Inv.A.add(sm, C, delta) +   (Inv.a.add(sm, C, delta) - Inv.A.add(sm, C, delta))/6,
                     Inv.A.add(sm, C, delta) +   (Inv.a.add(sm, C, delta) - Inv.A.add(sm, C, delta))/2,
                     Inv.A.add(sm, C, delta) + 5*(Inv.a.add(sm, C, delta) - Inv.A.add(sm, C, delta))/6,
                     Inv.a.add(sm, C, delta) + 0.1*(Inv.a.add(sm, C, delta) - Inv.A.add(sm, C, delta))
                    )
            }
        }
    if(h < 0.5){
        if(Inv.a.domRev(sm, h, C, delta) >= 1) {
        y  <-  c(Inv.A.domRev(sm, h, C, delta) - 0.1*(1 - Inv.A.domRev(sm, h, C, delta)),
                     Inv.A.domRev(sm, h, C, delta) +   (1 - Inv.A.domRev(sm, h, C, delta))/10,
                     Inv.A.domRev(sm, h, C, delta) +   (1 - Inv.A.domRev(sm, h, C, delta))/5,
                     Inv.A.domRev(sm, h, C, delta) +   (1 - Inv.A.domRev(sm, h, C, delta))/2,
                     1
                    )
            }
        if(Inv.a.domRev(sm, h, C, delta) < 1) {
            y  <-  c(Inv.A.domRev(sm, h, C, delta) - 0.1*(Inv.a.domRev(sm, h, C, delta) - Inv.A.domRev(sm, h, C, delta)),
                     Inv.A.domRev(sm, h, C, delta) +   (Inv.a.domRev(sm, h, C, delta) - Inv.A.domRev(sm, h, C, delta))/6,
                     Inv.A.domRev(sm, h, C, delta) +   (Inv.a.domRev(sm, h, C, delta) - Inv.A.domRev(sm, h, C, delta))/2,
                     Inv.A.domRev(sm, h, C, delta) + 5*(Inv.a.domRev(sm, h, C, delta) - Inv.A.domRev(sm, h, C, delta))/6,
                     Inv.a.domRev(sm, h, C, delta) + 0.1*(Inv.a.domRev(sm, h, C, delta) - Inv.A.domRev(sm, h, C, delta))
                    )
            }
        }
    y[y > 1]  <-  1
    y[y < 0]  <-  0
    y
}


##############################################################
##############################################################
##  Final figures for paper

#' Multipanel plot summarizing results for invasion of sterility alleles into 
#' populations at 1-locus EQ for the SA locus
#'
#' @title Figure 1
#' @author Colin Olito.
#' @export
Fig.1  <-  function() {

    # Import data
    data1  <-  read.csv("./output/data/EQInvAnalyses/Gyn-ObOut-Add-EQInv.csv", header=TRUE)
    data2  <-  read.csv("./output/data/EQInvAnalyses/Gyn-partSelf-C25-delta80-strgSel-Add-EQInv.csv", header=TRUE)
    data3  <-  read.csv("./output/data/EQInvAnalyses/Gyn-partSelf-C75-delta20-strgSel-Add-EQInv.csv", header=TRUE)
    data4  <-  read.csv("./output/data/EQInvAnalyses/And-ObOut-Add-EQInv.csv", header=TRUE)

    # Calculate Pr(Inv) for each data set
    d1  <-  ddply(data1, ~ k*r, summarize, pInv=(sum(DiffEQInvEig[DiffEQInvEig != 0])/length(DiffEQInvEig)))
    d2  <-  ddply(data2, ~ k*r, summarize, pInv=(sum(DiffEQInvEig[DiffEQInvEig != 0])/length(DiffEQInvEig)))
    d3  <-  ddply(data3, ~ k*r, summarize, pInv=(sum(DiffEQInvEig[DiffEQInvEig != 0])/length(DiffEQInvEig)))
    d4  <-  ddply(data4, ~ k*r, summarize, pInv=(sum(DiffEQInvEig[DiffEQInvEig != 0])/length(DiffEQInvEig)))
    # correct inconclusive Eigenvalue evaluations
    d1$pInv[d1$k == max(d1$k)]  <-  1
    d2$pInv[d2$k == max(d2$k)]  <-  1
    d4$pInv[d4$k == max(d4$k)]  <-  1

    # Color scheme
    COLS  <-  c(transparentColor('dodgerblue4', opacity=0.75),
                transparentColor('darkolivegreen', opacity=0.75),
                transparentColor('tomato2', opacity=0.75),
                transparentColor('dodgerblue2', opacity=0.75))
    COLS.bg  <-  c(transparentColor('dodgerblue4', opacity=0.5),
                transparentColor('darkolivegreen', opacity=0.5),
                transparentColor('tomato2', opacity=0.5),
                transparentColor('dodgerblue2', opacity=0.5))

    # k index for easy plotting
    ks1  <-  unique(d1$k)
    ks2  <-  unique(d2$k)
    ks3  <-  unique(d3$k)
    ks4  <-  unique(d4$k)

    # Set plot layout
    layout.mat <- matrix(c(1:6), nrow=3, ncol=2, byrow=FALSE)
    layout <- layout(layout.mat,respect=TRUE)

##  Column 1: Gynodioecy
    ##  Panel 1: Obligate Outcrossing
        # Calculate pretty x-values for plotting
        x  <-  unique(d1$r)
        xgap  <-  ifelse(x>0.2,x-0.3,x)
        xat  <-  pretty(xgap)
        xlab  <-  c(0.00, 0.01, 0.02, 0.10, 0.50)
        # make plot
        par(omi=rep(0.3, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(xgap)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot points
        lines(pInv[k==ks1[1]] ~ xat, col=COLS[3], lwd=2, data=d1)
        lines(pInv[k==ks1[2]] ~ xat, col=COLS[2], lwd=2, data=d1)
        lines(pInv[k==ks1[3]] ~ xat, col=COLS[1], lwd=2, data=d1)
        points(pInv[k==ks1[1]] ~ xat, pch=21, col=COLS[3], cex=1, bg=COLS.bg[3], data=d1)
        points(pInv[k==ks1[2]] ~ xat, pch=21, col=COLS[2], cex=1, bg=COLS.bg[2], data=d1)
        points(pInv[k==ks1[3]] ~ xat, pch=21, col=COLS[1], cex=1, bg=COLS.bg[1], data=d1)
        # axes
        axis(1, las=1, at=xat, labels=NA)
        axis(2, las=1)
        axis.break(1,0.175)
        axis.break(1,0.125)
        # Plot labels etc.
        proportionalLabel(0.5, 1.25, expression(paste("Gynodioecy")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(C," = 0")), cex=1, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.325, 0.5, expression(paste(Pr(inv))), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        legend(
              x       =  usr[2]*0.5,
              y       =  usr[4]*0.35,
              legend  =  c(
                          expression(paste(italic(k)~"="~1.0)),
                          expression(paste(italic(k)~"="~0.875)),
                          expression(paste(italic(k)~"="~0.75))),
              pch     =  c(21,21,21),
              pt.bg   =  c(COLS[1],COLS[2],COLS[3]),
              col     =  c(COLS[1],COLS[2],COLS[3]),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )

    ##  Panel 2: Low Selfing, High Inbreeding Depression
        # Calculate pretty x-values for plotting
        x  <-  unique(d2$r)
        xgap  <-  ifelse(x>0.2,x-0.3,x)
        xat  <-  pretty(xgap)
        xlab  <-  c(0.00, 0.01, 0.02, 0.10, 0.50)
        # make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(xgap)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot points
        lines(pInv[k==ks2[1]] ~ xat, col=COLS[4], lwd=2, data=d2)
        lines(pInv[k==ks2[2]] ~ xat, col=COLS[3], lwd=2, data=d2)
        lines(pInv[k==ks2[3]] ~ xat, col=COLS[2], lwd=2, data=d2)
        lines(pInv[k==ks2[4]] ~ xat, col=COLS[1], lwd=2, data=d2)
        points(pInv[k==ks2[1]] ~ xat, pch=21, col=COLS[4], cex=1, bg=COLS.bg[4], data=d2)
        points(pInv[k==ks2[2]] ~ xat, pch=21, col=COLS[3], cex=1, bg=COLS.bg[3], data=d2)
        points(pInv[k==ks2[3]] ~ xat, pch=21, col=COLS[2], cex=1, bg=COLS.bg[2], data=d2)
        points(pInv[k==ks2[4]] ~ xat, pch=21, col=COLS[1], cex=1, bg=COLS.bg[1], data=d2)
        # axes
        axis(1, las=1, at=xat, labels=NA)
        axis(2, las=1)
        axis.break(1,0.175)
        axis.break(1,0.125)
        # Plot labels etc.
        proportionalLabel(0.05, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(C," = 1/4, ",delta," = 4/5")), cex=1, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.325, 0.5, expression(paste(Pr(inv))), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
#        legend(
#              x       =  usr[2]*0.45,
#              y       =  usr[4]*0.35,
#              legend  =  c(
#                          expression(paste(hat(k))),
#                          expression(paste(hat(k)-0.1)),
#                          expression(paste(hat(k)-0.2)),
#                          expression(paste(hat(k)-0.3))),
#              pch     =  c(21,21,21),
#              pt.bg   =  c(COLS[1],COLS[2],COLS[3],COLS[4]),
#              col     =  c(COLS[1],COLS[2],COLS[3],COLS[4]),
#              cex     =  1,
#              xjust   =  1,
#              yjust   =  1,
#              bty     =  'n',
#              border  =  NA
#    )

    ##  Panel 3: High Selfing, Low Inbreeding Depression
        # Calculate pretty x-values for plotting
        x  <-  unique(d3$r)
        xgap  <-  ifelse(x>0.2,x-0.3,x)
        xat  <-  pretty(xgap)
        xlab  <-  c(0.00, 0.01, 0.02, 0.10, 0.50)
        # make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(xgap)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot points
        lines(pInv[k==ks3[1]] ~ xat, col=COLS[4], lwd=2, data=d3)
        lines(pInv[k==ks3[2]] ~ xat, col=COLS[3], lwd=2, data=d3)
        lines(pInv[k==ks3[3]] ~ xat, col=COLS[2], lwd=2, data=d3)
        lines(pInv[k==ks3[4]] ~ xat, col=COLS[1], lwd=2, data=d3)
        points(pInv[k==ks3[1]] ~ xat, pch=21, col=COLS[4], cex=1, bg=COLS.bg[4], data=d3)
        points(pInv[k==ks3[2]] ~ xat, pch=21, col=COLS[3], cex=1, bg=COLS.bg[3], data=d3)
        points(pInv[k==ks3[3]] ~ xat, pch=21, col=COLS[2], cex=1, bg=COLS.bg[2], data=d3)
        points(pInv[k==ks3[4]] ~ xat, pch=21, col=COLS[1], cex=1, bg=COLS.bg[1], data=d3)
        # axes
        axis(1, las=1, at=xat, labels=xlab)
        axis(2, las=1)
        axis.break(1,0.175)
        axis.break(1,0.125)
        # Plot labels etc.
        proportionalLabel(0.05, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  1.075, expression(paste(C," = 3/4, ",delta," = 1/5")), cex=1, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.325, 0.5, expression(paste(Pr(inv))), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.3, expression(paste(italic(r))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        legend(
              x       =  usr[2]*0.975,
              y       =  usr[4]*1.03,
              legend  =  c(
                          expression(paste(hat(italic(k)))),
                          expression(paste(hat(italic(k))-0.1)),
                          expression(paste(hat(italic(k))-0.2)),
                          expression(paste(hat(italic(k))-0.3))),
              pch     =  c(21,21,21),
              pt.bg   =  c(COLS[1],COLS[2],COLS[3],COLS[4]),
              col     =  c(COLS[1],COLS[2],COLS[3],COLS[4]),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )

##  Column 2: Androdioecy
    ##  Panel 4: Obligate Outcrossing
        # Calculate pretty x-values for plotting
        x  <-  unique(d4$r)
        xgap  <-  ifelse(x>0.2,x-0.3,x)
        xat  <-  pretty(xgap)
        xlab  <-  c(0.00, 0.01, 0.02, 0.10, 0.50)
        # make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(xgap)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot points
        lines(pInv[k==ks1[1]] ~ xat, col=COLS[3], lwd=2, data=d4)
        lines(pInv[k==ks1[2]] ~ xat, col=COLS[2], lwd=2, data=d4)
        lines(pInv[k==ks1[3]] ~ xat, col=COLS[1], lwd=2, data=d4)
        points(pInv[k==ks1[1]] ~ xat, pch=21, col=COLS[3], cex=1, bg=COLS.bg[3], data=d4)
        points(pInv[k==ks1[2]] ~ xat, pch=21, col=COLS[2], cex=1, bg=COLS.bg[2], data=d4)
        points(pInv[k==ks1[3]] ~ xat, pch=21, col=COLS[1], cex=1, bg=COLS.bg[1], data=d4)
        # axes
        axis(1, las=1, at=xat, labels=NA)
        axis(2, las=1, labels=NA)
        axis.break(1,0.175)
        axis.break(1,0.125)
        # Plot labels etc.
        proportionalLabel(0.5, 1.25, expression(paste("Androdioecy")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(C," = 0")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    ##  Panel 5: Low Selfing, High Inbreeding Depression
        # Calculate pretty x-values for plotting
        x  <-  unique(d5$r)
        xgap  <-  ifelse(x>0.2,x-0.3,x)
        xat  <-  pretty(xgap)
        xlab  <-  c(0.00, 0.01, 0.02, 0.10, 0.50)
        # make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(xgap)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot points
        lines(pInv[k==ks5[1]] ~ xat, col=COLS[4], lwd=2, data=d5)
        lines(pInv[k==ks5[2]] ~ xat, col=COLS[3], lwd=2, data=d5)
        lines(pInv[k==ks5[3]] ~ xat, col=COLS[2], lwd=2, data=d5)
        lines(pInv[k==ks5[4]] ~ xat, col=COLS[1], lwd=2, data=d5)
        points(pInv[k==ks5[1]] ~ xat, pch=21, col=COLS[4], cex=1, bg=COLS.bg[4], data=d5)
        points(pInv[k==ks5[2]] ~ xat, pch=21, col=COLS[3], cex=1, bg=COLS.bg[3], data=d5)
        points(pInv[k==ks5[3]] ~ xat, pch=21, col=COLS[2], cex=1, bg=COLS.bg[2], data=d5)
        points(pInv[k==ks5[4]] ~ xat, pch=21, col=COLS[1], cex=1, bg=COLS.bg[1], data=d5)
        # axes
        axis(1, las=1, at=xat, labels=NA)
        axis(2, las=1, labels=NA)
        axis.break(1,0.175)
        axis.break(1,0.125)
        # Plot labels etc.
        proportionalLabel(0.05, 1.075, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(C," = 1/4, ",delta," = 0.8")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        legend(
#              x       =  usr[2]*0.45,
#              y       =  usr[4]*0.35,
#              legend  =  c(
#                          expression(paste(hat(k))),
#                          expression(paste(hat(k)-0.1)),
#                          expression(paste(hat(k)-0.2)),
#                          expression(paste(hat(k)-0.3))),
#              pch     =  c(21,21,21),
#              pt.bg   =  c(COLS[1],COLS[2],COLS[3],COLS[4]),
#              col     =  c(COLS[1],COLS[2],COLS[3],COLS[4]),
#              cex     =  1,
#              xjust   =  1,
#              yjust   =  1,
#              bty     =  'n',
#              border  =  NA
#    )

    ##  Panel 6: High Selfing, Low Inbreeding Depression
        # Calculate pretty x-values for plotting
        x  <-  unique(d6$r)
        xgap  <-  ifelse(x>0.2,x-0.3,x)
        xat  <-  pretty(xgap)
        xlab  <-  c(0.00, 0.01, 0.02, 0.10, 0.50)
        # make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(xgap)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot points
        lines(pInv[k==ks6[1]] ~ xat, col=COLS[4], lwd=2, data=d6)
        lines(pInv[k==ks6[2]] ~ xat, col=COLS[3], lwd=2, data=d6)
        lines(pInv[k==ks6[3]] ~ xat, col=COLS[2], lwd=2, data=d6)
        lines(pInv[k==ks6[4]] ~ xat, col=COLS[1], lwd=2, data=d6)
        points(pInv[k==ks6[1]] ~ xat, pch=21, col=COLS[4], cex=1, bg=COLS.bg[4], data=d6)
        points(pInv[k==ks6[2]] ~ xat, pch=21, col=COLS[3], cex=1, bg=COLS.bg[3], data=d6)
        points(pInv[k==ks6[3]] ~ xat, pch=21, col=COLS[2], cex=1, bg=COLS.bg[2], data=d6)
        points(pInv[k==ks6[4]] ~ xat, pch=21, col=COLS[1], cex=1, bg=COLS.bg[1], data=d6)
        # axes
        axis(1, las=1, at=xat, labels=xlab)
        axis(2, las=1, labels=NA)
        axis.break(1,0.175)
        axis.break(1,0.125)
        # Plot labels etc.
        proportionalLabel(0.05, 1.075, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(C," = 1/4, ",delta," = 0.2")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.35, expression(paste(italic(r))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        legend(
#              x       =  usr[2]*0.45,
#              y       =  usr[4]*0.35,
#              legend  =  c(
#                          expression(paste(hat(k))),
#                          expression(paste(hat(k)-0.1)),
#                          expression(paste(hat(k)-0.2)),
#                          expression(paste(hat(k)-0.3))),
#              pch     =  c(21,21,21),
#              pt.bg   =  c(COLS[1],COLS[2],COLS[3],COLS[4]),
#              col     =  c(COLS[1],COLS[2],COLS[3],COLS[4]),
#              cex     =  1,
#              xjust   =  1,
#              yjust   =  1,
#              bty     =  'n',
#              border  =  NA
#    )

}



##############################################################
##############################################################
##  Exploratory plots

#' Exploratory plots: Invasion of dominant sterility allele into populations at 1-locus SA equilibrium. 
#'                    -- Obligate outcrossing
#'                    -- Additive fitness effects at SA locus
#' @title Invasion of dominant sterility allele into populations
#' @author Colin Olito
#' @export
EQInv.ObOut.Add  <-  function(df="./output/data/EQInvAnalyses/Gyn-ObOut-Add-EQInv.csv") {

    # Import data
    data  <-  read.csv(df, header=TRUE)

    # Color scheme
    COLS  <-  transparentColor('dodgerblue4', opacity=0.2)

    # Calculate 1-locus SA invasion criteria to illustrate 
    # boundaries for polymorphic populations
    sms  <-  seq(0,1,by=0.0001)
    Ainv  <-  Inv.A.add(sms, C=0, delta=0)
    Ainv[Ainv > 1]  <-  1.00001
    ainv  <-  Inv.a.add(sms, C=0, delta=0)

    # k index for easy plotting
    ks  <-  unique(data$k)

    # Set plot layout
    layout.mat <- matrix(c(1:15), nrow=5, ncol=3, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

##  Row 1: r = 0
    ##  Panel 1: k = 1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==0 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==0]), precision=3)
        # Make plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==1 & r==0 & DiffEQInvEig != 0] ~ sm[k==1 & r==0 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.5, 1.25, expression(paste(italic(k), ' = ', 1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 2: k = 0.875
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==0 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==0]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==0 & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==0 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.5, 1.25, expression(paste(italic(k), ' = ', 0.875)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 3: k = 0.75
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==0 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==0]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==0 & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==0 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.5, 1.25, expression(paste(italic(k), ' = ', 0.75)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)



##  Row 2: r = 0.01
    ##  Panel 4: k = 1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.01 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.01]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==1 & r==0.01 & DiffEQInvEig != 0] ~ sm[k==1 & r==0.01 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.01")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 5: k = 0.875
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.01 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.01]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==0.01 & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==0.01 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.05, 1.075, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 6: k = 0.75
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.01 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.01]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==0.01 & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==0.01 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.05, 1.075, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)


##  Row 3: r = 0.02
    ##  Panel 7: k = 1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.02 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.02]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==1 & r==0.02 & DiffEQInvEig != 0] ~ sm[k==1 & r==0.02 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.02")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 8: k = 0.875
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.02 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.02]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==0.02 & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==0.02 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.05, 1.075, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 9: k = 0.75
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.02 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.02]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==0.02 & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==0.02 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.05, 1.075, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)


##  Row 4: r = 0.1
    ##  Panel 10: k = 1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.1 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.1]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==1 & r==0.1 & DiffEQInvEig != 0] ~ sm[k==1 & r==0.1 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.1")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'J', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 11: k = 0.875
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.1 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.1]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==0.1 & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==0.1 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.05, 1.075, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 12: k = 0.75
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.1 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.1]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==0.1 & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==0.1 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.05, 1.075, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)


##  Row 5: r = 0.5
    ##  Panel 13: k = 1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.5 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.5]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==1 & r==0.5 & DiffEQInvEig != 0] ~ sm[k==1 & r==0.5 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.5")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'M', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 14: k = 0.875
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.5 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.5]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==0.5 & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==0.5 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'N', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 15: k = 0.75
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.5 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.5]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==0.5 & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==0.5 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'O', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)
}




#' Exploratory plots: Invasion of dominant sterility allele into populations at 1-locus SA equilibrium. 
#'                    -- Partial Selfing
#'                    -- Additive fitness effects at SA locus
#' @title Invasion of dominant male sterility allele into populations
#' @author Colin Olito
#' @export
EQInv.PartSelf.Add  <-  function(df="./output/data/EQInvAnalyses/Gyn-partSelf-C25-delta20-strgSel-Add-EQInv.csv") {

    # Import data
    data  <-  read.csv(df, header=TRUE)

    # Color scheme
    COLS  <-  transparentColor('dodgerblue4', opacity=0.2)

    # Calculate 1-locus SA invasion criteria to illustrate 
    # boundaries for polymorphic populations
    sms  <-  seq(0,1,by=0.01)
    Ainv  <-  Inv.A.add(sms, C = data$C[1], delta = data$delta[1])
    Ainv[Ainv > 1]  <-  1.00001
    ainv  <-  Inv.a.add(sms, C = data$C[1], delta = data$delta[1])

    # k index for easy plotting
    ks  <-  unique(data$k)

    # Set plot layout
    layout.mat <- matrix(c(1:20), nrow=5, ncol=4, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

##  Row 1: r = 0
    ##  Panel 1: k = kCrit
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==0 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==0]), precision=3)
        # Make plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[1] & r==0 & DiffEQInvEig != 0] ~ sm[k==ks[1] & r==0 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.5, 1.25, expression(paste(hat(italic(k))," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.7, 1.23, substitute(k,list(k=ks[1])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 2: k = kCrit - 0.1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==0 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==0]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==0 & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==0 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(hat(italic(k)) - 0.1," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.88, 1.23, substitute(k,list(k=ks[2])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 3: k = kCrit - 0.2
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==0 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==0]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==0 & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==0 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(hat(italic(k)) - 0.2," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.88, 1.23, substitute(k,list(k=ks[3])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 4: k = kCrit - 0.3
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[4] & data$r==0 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[4] & data$r==0]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[4] & r==0 & DiffEQInvEig != 0] ~ sm[k==ks[4] & r==0 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(hat(italic(k)) - 0.3," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.88, 1.23, substitute(k,list(k=ks[4])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)


##  Row 2: r = 0.01
    ##  Panel 5: k = kCrit
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.01 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.01]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[1] & r==0.01 & DiffEQInvEig != 0] ~ sm[k==ks[1] & r==0.01 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.01")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 6: k = kCrit - 0.1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.01 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.01]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==0.01 & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==0.01 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 7: k = kCrit - 0.2
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.01 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.01]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==0.01 & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==0.01 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 8: k = kCrit - 0.3
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[4] & data$r==0.01 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[4] & data$r==0.01]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[4] & r==0.01 & DiffEQInvEig != 0] ~ sm[k==ks[4] & r==0.01 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

##  Row 3: r = 0.02
    ##  Panel 9: k = kCrit
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.02 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.02]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[1] & r==0.02 & DiffEQInvEig != 0] ~ sm[k==ks[1] & r==0.02 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.02")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 10: k = kCrit - 0.1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.02 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.02]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==0.02 & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==0.02 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 11: k = kCrit - 0.2
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.02 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.02]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==0.02 & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==0.02 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'J', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 12: k = kCrit - 0.3
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[4] & data$r==0.02 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[4] & data$r==0.02]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[4] & r==0.02 & DiffEQInvEig != 0] ~ sm[k==ks[4] & r==0.02 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)


##  Row 4: r = 0.1
    ##  Panel 13: k = kCrit
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.1 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.1]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[1] & r==0.1 & DiffEQInvEig != 0] ~ sm[k==ks[1] & r==0.1 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.1")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 14: k = kCrit - 0.1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.1 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.1]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==0.1 & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==0.1 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'M', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 15: k = kCrit - 0.2
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.1 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.1]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==0.1 & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==0.1 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'N', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 16: k = kCrit - 0.3
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[4] & data$r==0.1 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[4] & data$r==0.1]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[4] & r==0.1 & DiffEQInvEig != 0] ~ sm[k==ks[4] & r==0.1 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'O', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)


##  Row 5: r = 0.5
    ##  Panel 17: k = kCrit
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.5 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==0.5]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[1] & r==0.5 & DiffEQInvEig != 0] ~ sm[k==ks[1] & r==0.5 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.5")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'P', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 18: k = kCrit - 0.1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.5 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==0.5]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==0.5 & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==0.5 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'Q', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 19: k = kCrit - 0.2
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.5 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==0.5]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==0.5 & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==0.5 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'R', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 20: k = kCrit - 0.3
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[4] & data$r==0.5 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[4] & data$r==0.5]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[4] & r==0.5 & DiffEQInvEig != 0] ~ sm[k==ks[4] & r==0.5 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'S', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

}