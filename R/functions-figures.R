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
    data1  <-  read.csv("./output/data/EQInvAnalyses/Gyn-wksel-ObOut-Add-EQInv.csv", header=TRUE)
    data2  <-  read.csv("./output/data/EQInvAnalyses/Gyn-wksel-partSelf-C25-delta80-Add-EQInv.csv", header=TRUE)
    data3  <-  read.csv("./output/data/EQInvAnalyses/Gyn-wksel-partSelf-C75-delta20-Add-EQInv.csv", header=TRUE)
    data4  <-  read.csv("./output/data/EQInvAnalyses/And-wksel-ObOut-Add-EQInv.csv", header=TRUE)
    data5  <-  read.csv("./output/data/EQInvAnalyses/And-wksel-partSelf-C25-delta80-Add-EQInv.csv", header=TRUE)
    data6  <-  read.csv("./output/data/EQInvAnalyses/And-wksel2-partSelf-C75-delta20-Add-EQInv.csv", header=TRUE)

    # k index for easy plotting
    ks1  <-  unique(data1$k) 
    ks2  <-  unique(data2$k)
    ks3  <-  unique(data3$k)
    ks4  <-  unique(data4$k)
    ks5  <-  unique(data5$k)
    ks6  <-  unique(data6$k)

    # Calculate Pr(Inv) for each data set
    d1  <-  ddply(data1, ~ k*r, summarize, pInv=(length(DiffEQInvEig[DiffEQInvEig != 0])/length(DiffEQInvEig)))
    d2  <-  ddply(data2, ~ k*r, summarize, pInv=(length(DiffEQInvEig[DiffEQInvEig != 0])/length(DiffEQInvEig)))
    d3  <-  ddply(data3, ~ k*r, summarize, pInv=(length(DiffEQInvEig[DiffEQInvEig != 0])/length(DiffEQInvEig)))
    d4  <-  ddply(data4, ~ k*r, summarize, pInv=(length(DiffEQInvEig[DiffEQInvEig != 0])/length(DiffEQInvEig)))
    d5  <-  ddply(data5, ~ k*r, summarize, pInv=(length(DiffEQInvEig[DiffEQInvEig != 0])/length(DiffEQInvEig)))
    d6  <-  ddply(data6, ~ k*r, summarize, pInv=(length(DiffEQInvEig[DiffEQInvEig != 0])/length(DiffEQInvEig)))

    # Color scheme
    COLS  <-  c(transparentColor('dodgerblue4', opacity=0.75),
                transparentColor('darkolivegreen', opacity=0.75),
                transparentColor('tomato2', opacity=0.75),
                transparentColor('dodgerblue2', opacity=0.75))
    COLS.bg  <-  c(transparentColor('dodgerblue4', opacity=0.5),
                transparentColor('darkolivegreen', opacity=0.5),
                transparentColor('tomato2', opacity=0.5),
                transparentColor('dodgerblue2', opacity=0.5))

    # Set plot layout
    layout.mat <- matrix(c(1:6), nrow=3, ncol=2, byrow=FALSE)
    layout <- layout(layout.mat,respect=TRUE)

##  Column 1: Gynodioecy
    ##  Panel 1: Obligate Outcrossing
        # Calculate pretty x-values for plotting
        x  <-  unique(d1$r)
        xat  <-  pretty(x)
        xlab  <-  c(0.00, 0.01, 0.05, 0.10, 0.20, 0.50)
        # make plot
        par(omi=rep(0.3, 4), mar = c(3.5,3.5,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(x)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot points
        lines(pInv[k==ks1[1]]  ~ xat, col=COLS[4], lwd=2, data=d1)
        lines(pInv[k==ks1[2]]  ~ xat, col=COLS[3], lwd=2, data=d1)
        lines(pInv[k==ks1[3]]  ~ xat, col=COLS[2], lwd=2, data=d1)
        lines(pInv[k==ks1[4]]  ~ xat, col=COLS[1], lwd=2, data=d1)
        points(pInv[k==ks1[1]] ~ xat, pch=21, col=COLS[4], cex=1, bg=COLS.bg[4], data=d1)
        points(pInv[k==ks1[2]] ~ xat, pch=21, col=COLS[3], cex=1, bg=COLS.bg[3], data=d1)
        points(pInv[k==ks1[3]] ~ xat, pch=21, col=COLS[2], cex=1, bg=COLS.bg[2], data=d1)
        points(pInv[k==ks1[4]] ~ xat, pch=21, col=COLS[1], cex=1, bg=COLS.bg[1], data=d1)
        # axes
        axis(1, las=1, at=xat, labels=NA)
        axis(2, las=1)
        axis.break(1,0.15)
        axis.break(1,0.25)
        axis.break(1,0.35)
        axis.break(1,0.45)
        # Plot labels etc.
        proportionalLabel(0.5, 1.25, expression(paste("Gynodioecy")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(C," = 0")), cex=1, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.325, 0.5, expression(paste(Fraction~of~parameter~space)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        #legend
        legend(
              x       =  usr[2]*0.48,
              y       =  usr[4]*0.39,
              legend  =  c(
                          expression(paste(italic(k)~"="~italic(hat(k))%*%1.1)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.99)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.95)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.9))),
              pch     =  c(21,21,21,21),
              pt.bg   =  c(COLS[4],COLS[3],COLS[2],COLS[1]),
              col     =  c(COLS[4],COLS[3],COLS[2],COLS[1]),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )

    ##  Panel 2: Low Selfing, High Inbreeding Depression
        # Calculate pretty x-values for plotting
        x     <-  unique(d2$r)
        xat   <-  pretty(x)
        xlab  <-  c(0.00, 0.01, 0.02, 0.10, 0.50)
        # make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(x)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
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
        axis(1, las=1, at=xat, labels=x)
        axis(2, las=1)
        axis.break(1,0.15)
        axis.break(1,0.25)
        axis.break(1,0.35)
        axis.break(1,0.45)
        # Plot labels etc.
        proportionalLabel(0.05, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(C," = ",0.25,", ",delta," = 0.8")), cex=1, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.325, 0.5, expression(paste(Fraction~of~parameter~space)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)

    ##  Panel 3: High Selfing, Low Inbreeding Depression
        # Calculate pretty x-values for plotting
        x     <-  unique(d3$r)
        xat   <-  pretty(x)
        xlab  <-  c(0.00, 0.01, 0.02, 0.10, 0.50)
        # make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(x)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
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
        axis(1, las=1, at=xat, labels=x)
        axis(2, las=1)
        axis.break(1,0.15)
        axis.break(1,0.25)
        axis.break(1,0.35)
        axis.break(1,0.45)
        # Plot labels etc.
        proportionalLabel(0.05, 1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  1.075, expression(paste(C," = ",0.75,", ",delta, " = ",0.2)), cex=1, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.325, 0.5, expression(paste(Fraction~of~parameter~space)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.3, expression(paste(italic(r))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)

##  Column 2: Androdioecy
    ##  Panel 4: Obligate Outcrossing
        # Calculate pretty x-values for plotting
        x  <-  unique(d4$r)
        xat  <-  pretty(x)
        xlab  <-  c(0.00, 0.01, 0.05, 0.10, 0.20, 0.50)
        # make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(x)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot points
        lines(pInv[k==ks1[1]]  ~ xat, col=COLS[4], lwd=2, data=d4)
        lines(pInv[k==ks1[2]]  ~ xat, col=COLS[3], lwd=2, data=d4)
        lines(pInv[k==ks1[3]]  ~ xat, col=COLS[2], lwd=2, data=d4)
        lines(pInv[k==ks1[4]]  ~ xat, col=COLS[1], lwd=2, data=d4)
        points(pInv[k==ks1[1]] ~ xat, pch=21, col=COLS[4], cex=1, bg=COLS.bg[4], data=d4)
        points(pInv[k==ks1[2]] ~ xat, pch=21, col=COLS[3], cex=1, bg=COLS.bg[3], data=d4)
        points(pInv[k==ks1[3]] ~ xat, pch=21, col=COLS[2], cex=1, bg=COLS.bg[2], data=d4)
        points(pInv[k==ks1[4]] ~ xat, pch=21, col=COLS[1], cex=1, bg=COLS.bg[1], data=d4)
        # axes
        axis(1, las=1, at=xat, labels=NA)
        axis(2, las=1, labels=NA)
        axis.break(1,0.15)
        axis.break(1,0.25)
        axis.break(1,0.35)
        axis.break(1,0.45)
        # Plot labels etc.
        proportionalLabel(0.5, 1.25, expression(paste("Androdioecy")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(C," = 0")), cex=1, adj=c(0.5, 0.5), xpd=NA)

    ##  Panel 5: Low Selfing, High Inbreeding Depression
        # Calculate pretty x-values for plotting
        x  <-  unique(d5$r)
        xat  <-  pretty(x)
        xlab  <-  c(0.00, 0.01, 0.02, 0.10, 0.50)
        # make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(x)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot points
        lines(pInv[k==ks5[1]]  ~ xat, col=COLS[4], lwd=2, data=d5)
        lines(pInv[k==ks5[2]]  ~ xat, col=COLS[3], lwd=2, data=d5)
        lines(pInv[k==ks5[3]]  ~ xat, col=COLS[2], lwd=2, data=d5)
        lines(pInv[k==ks5[4]]  ~ xat, col=COLS[1], lwd=2, data=d5)
        points(pInv[k==ks5[1]] ~ xat, pch=21, col=COLS[4], cex=1, bg=COLS.bg[4], data=d5)
        points(pInv[k==ks5[2]] ~ xat, pch=21, col=COLS[3], cex=1, bg=COLS.bg[3], data=d5)
        points(pInv[k==ks5[3]] ~ xat, pch=21, col=COLS[2], cex=1, bg=COLS.bg[2], data=d5)
        points(pInv[k==ks5[4]] ~ xat, pch=21, col=COLS[1], cex=1, bg=COLS.bg[1], data=d5)
        # axes
        axis(1, las=1, at=xat, labels=NA)
        axis(2, las=1, labels=NA)
        axis.break(1,0.15)
        axis.break(1,0.25)
        axis.break(1,0.35)
        axis.break(1,0.45)
        # Plot labels etc.
        proportionalLabel(0.05, 1.075, expression(paste(bold(E))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(C," = ",0.25,", ",delta," = 0.8")), cex=1, adj=c(0.5, 0.5), xpd=NA)

    ##  Panel 6: High Selfing, Low Inbreeding Depression
        x     <-  unique(d6$r)
        xat   <-  pretty(x)
        # make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(x)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
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
        axis(1, las=1, at=xat, labels=x)
        axis(2, las=1, labels=NA)
        axis.break(1,0.15)
        axis.break(1,0.25)
        axis.break(1,0.35)
        axis.break(1,0.45)
#        # Plot labels etc.
        proportionalLabel(0.05, 1.075, expression(paste(bold(F))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  1.075, expression(paste(C," = ",0.75,", ",delta, " = ",0.2)), cex=1, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.35, expression(paste(italic(r))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

}



#' Multipanel plot comparing 2-locus recursion simulation results
#' with 1-locus predictions
#'
#' @title Figure 2 alternate
#' @author Colin Olito.
#' @export
Fig2Alt  <-  function(dfGyn, dfAnd) {
    # Import data
    dataG  <-  read.csv(dfGyn, header=TRUE)
    dataA  <-  read.csv(dfAnd, header=TRUE)

    # Calculate equilibrium frequencies of M2, females, A, a
    dataG$q.m2        <-  (dataG$'F.12'/2) + (dataG$'F.14'/2) + dataG$'F.22' + (dataG$'F.23'/2) + dataG$'F.24' + (dataG$'F.34'/2) + dataG$'F.44' + 
                         (dataG$'G.12'/2) + (dataG$'G.14'/2) + dataG$'G.22' + (dataG$'G.23'/2) + dataG$'G.24' + (dataG$'G.34'/2) + dataG$'G.44'
    dataG$females     <-  dataG$'F.22' + dataG$'F.24' + dataG$'F.44' + dataG$'G.22' + dataG$'G.24' + dataG$'G.44'
    dataG$p.A         <-  dataG$'F.11' + dataG$'F.12' + (dataG$'F.13'/2) + (dataG$'F.14'/2) + dataG$'F.22' + (dataG$'F.23'/2) + (dataG$'F.24'/2) + 
                         dataG$'G.11' + dataG$'G.12' + (dataG$'G.13'/2) + (dataG$'G.14'/2) + dataG$'G.22' + (dataG$'G.23'/2) + (dataG$'G.24'/2)
    dataG$q.a         <-  (dataG$'F.13'/2) + (dataG$'F.14'/2) + (dataG$'F.23'/2) + (dataG$'F.24'/2) + dataG$'F.33' + dataG$'F.34' + dataG$'F.44' + 
                         (dataG$'G.13'/2) + (dataG$'G.14'/2) + (dataG$'G.23'/2) + (dataG$'G.24'/2) + dataG$'G.33' + dataG$'G.34' + dataG$'G.44' 
    dataG$diffFemales  <-  (dataG$females - dataG$ZHat)

    # Calculate equilibrium frequencies of M2, males, A, a
    dataA$q.m2        <-  (dataA$'F.12'/2) + (dataA$'F.14'/2) + dataA$'F.22' + (dataA$'F.23'/2) + dataA$'F.24' + (dataA$'F.34'/2) + dataA$'F.44' + 
                         (dataA$'G.12'/2) + (dataA$'G.14'/2) + dataA$'G.22' + (dataA$'G.23'/2) + dataA$'G.24' + (dataA$'G.34'/2) + dataA$'G.44'
    dataA$males     <-  dataA$'F.22' + dataA$'F.24' + dataA$'F.44' + dataA$'G.22' + dataA$'G.24' + dataA$'G.44'
    dataA$p.A         <-  dataA$'F.11' + dataA$'F.12' + (dataA$'F.13'/2) + (dataA$'F.14'/2) + dataA$'F.22' + (dataA$'F.23'/2) + (dataA$'F.24'/2) + 
                         dataA$'G.11' + dataA$'G.12' + (dataA$'G.13'/2) + (dataA$'G.14'/2) + dataA$'G.22' + (dataA$'G.23'/2) + (dataA$'G.24'/2)
    dataA$q.a         <-  (dataA$'F.13'/2) + (dataA$'F.14'/2) + (dataA$'F.23'/2) + (dataA$'F.24'/2) + dataA$'F.33' + dataA$'F.34' + dataA$'F.44' + 
                         (dataA$'G.13'/2) + (dataA$'G.14'/2) + (dataA$'G.23'/2) + (dataA$'G.24'/2) + dataA$'G.33' + dataA$'G.34' + dataA$'G.44' 
    dataA$diffMales  <-  (dataA$males - dataA$ZHat)

    # Color scheme
    COLS  <-  c(transparentColor('dodgerblue', opacity=0.8),
                transparentColor('dodgerblue4', opacity=0.8),
                transparentColor('darkolivegreen', opacity=0.8),
                transparentColor('tomato3', opacity=0.8),
                transparentColor('tomato', opacity=0.8))

    # k index for easy plotting
    rs  <-  unique(dataG$r)
    ks  <-  unique(dataG$k)

    # Set plot layout
    layout.mat <- matrix(c(1:2), nrow=1, ncol=2, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

### Equilibrium Female Frequencies
##  Panel 1: 
    dat  <-  subset(dataG, dataG$k == ks[1])
        par(omi=rep(0.3, 4), mar = c(3.5,3.5,2,2), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(0,max(dataG$females)*1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different 
        lines(ZHat[r==rs[1]] ~ C[r==rs[1]], lwd=2, col='black', cex=1, data=dat)
        lines(females[r==rs[1]] ~ C[r==rs[1]], lwd=2, col='grey70', lty=2, cex=1, data=dat)
        lines(females[r==rs[2]] ~ C[r==rs[2]], lwd=2, col='grey60', lty=2, cex=1, data=dat)
        lines(females[r==rs[3]] ~ C[r==rs[3]], lwd=2, col='grey50', lty=2, cex=1, data=dat)
        lines(females[r==rs[4]] ~ C[r==rs[4]], lwd=2, col='grey40', lty=2, cex=1, data=dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.5, 1.25, expression(paste("Gynodioecy")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(C))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35, 0.5, expression(paste("Eq. female frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        #legend
        legend(
              x       =  usr[2]*0.99,
              y       =  usr[4],
              legend  =  c(
                          expression(paste(italic(k)~"="~italic(hat(k))%*%1.1)),
                          expression(paste("1 locus"))),
              lty     =  c(2,1),
              lwd     =  2,
              col     =  'black',
              cex     =  0.75,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
    proportionalLabel(0.115, 0.92, expression(paste(italic(r)~"="~0.0)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.165, 0.709, expression(paste(italic(r)~"="~0.005)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.15, 0.6, expression(paste(italic(r)~"="~0.01)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.15, 0.482, expression(paste(italic(r)~"="~0.05)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    rm(dat)

### Equilibrium Male Frequencies
##  Panel 2: 
    dat  <-  subset(dataA, dataA$k == ks[1])
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(0,max(dataA$males)*1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different 
        lines(ZHat[r==rs[1]] ~ C[r==rs[1]], lwd=2, col='black', cex=1, data=dat)
        lines(males[r==rs[1]] ~ C[r==rs[1]], lwd=2, col='grey70', lty=2, cex=1, data=dat)
        lines(males[r==rs[2]] ~ C[r==rs[2]], lwd=2, col='grey60', lty=2, cex=1, data=dat)
        lines(males[r==rs[3]] ~ C[r==rs[3]], lwd=2, col='grey50', lty=2, cex=1, data=dat)
        lines(males[r==rs[4]] ~ C[r==rs[4]], lwd=2, col='grey40', lty=2, cex=1, data=dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.5, 1.25, expression(paste("Androdioecy")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35, 0.5, expression(paste("Eq. male frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.3, expression(paste(italic(C))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.115, 0.92, expression(paste(italic(r)~"="~0.0)),   cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.165, 0.705, expression(paste(italic(r)~"="~0.005)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.15, 0.618, expression(paste(italic(r)~"="~0.01)),  cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.15, 0.505, expression(paste(italic(r)~"="~0.05)),  cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    rm(dat)
}



#' Multipanel plot summarizing results recursion simulations
#'
#' @title Figure 2
#' @author Colin Olito.
#' @export
Fig.2  <-  function() {

    # Import data
    data   <-  read.csv(file="./output/data/simResults/gyn-recess_dStar0.8_a1_sm0.1_add.csv", header=TRUE)
    data2  <-  read.csv(file="./output/data/simResults/and-recess_dStar0.8_a1_sm0.1_add.csv", header=TRUE)

    # Calculate equilibrium frequencies of M2, females, A, a
    data$q.m2        <-  (data$'F.12'/2) + (data$'F.14'/2) + data$'F.22' + (data$'F.23'/2) + data$'F.24' + (data$'F.34'/2) + data$'F.44' + 
                         (data$'G.12'/2) + (data$'G.14'/2) + data$'G.22' + (data$'G.23'/2) + data$'G.24' + (data$'G.34'/2) + data$'G.44'
    data$females     <-  data$'F.22' + data$'F.24' + data$'F.44' + data$'G.22' + data$'G.24' + data$'G.44'
    data$p.A         <-  data$'F.11' + data$'F.12' + (data$'F.13'/2) + (data$'F.14'/2) + data$'F.22' + (data$'F.23'/2) + (data$'F.24'/2) + 
                         data$'G.11' + data$'G.12' + (data$'G.13'/2) + (data$'G.14'/2) + data$'G.22' + (data$'G.23'/2) + (data$'G.24'/2)
    data$q.a         <-  (data$'F.13'/2) + (data$'F.14'/2) + (data$'F.23'/2) + (data$'F.24'/2) + data$'F.33' + data$'F.34' + data$'F.44' + 
                         (data$'G.13'/2) + (data$'G.14'/2) + (data$'G.23'/2) + (data$'G.24'/2) + data$'G.33' + data$'G.34' + data$'G.44' 
    data$diffFemales  <-  (data$females - data$ZHat)

   # Calculate equilibrium frequencies of M2, males, A, a
    data2$q.m2        <-  (data2$'F.12'/2) + (data2$'F.14'/2) + data2$'F.22' + (data2$'F.23'/2) + data2$'F.24' + (data2$'F.34'/2) + data2$'F.44' + 
                         (data2$'G.12'/2) + (data2$'G.14'/2) + data2$'G.22' + (data2$'G.23'/2) + data2$'G.24' + (data2$'G.34'/2) + data2$'G.44'
    data2$Males       <-  data2$'F.22' + data2$'F.24' + data2$'F.44' + data2$'G.22' + data2$'G.24' + data2$'G.44'
    data2$p.A         <-  data2$'F.11' + data2$'F.12' + (data2$'F.13'/2) + (data2$'F.14'/2) + data2$'F.22' + (data2$'F.23'/2) + (data2$'F.24'/2) + 
                         data2$'G.11' + data2$'G.12' + (data2$'G.13'/2) + (data2$'G.14'/2) + data2$'G.22' + (data2$'G.23'/2) + (data2$'G.24'/2)
    data2$q.a         <-  (data2$'F.13'/2) + (data2$'F.14'/2) + (data2$'F.23'/2) + (data2$'F.24'/2) + data2$'F.33' + data2$'F.34' + data2$'F.44' + 
                         (data2$'G.13'/2) + (data2$'G.14'/2) + (data2$'G.23'/2) + (data2$'G.24'/2) + data2$'G.33' + data2$'G.34' + data2$'G.44' 
    data2$diffMales   <-  (data2$Males - data2$ZHat)

    # Color scheme
     COLS  <-  c(transparentColor('dodgerblue', opacity=0.8),
                transparentColor('dodgerblue4', opacity=0.8),
                transparentColor('darkolivegreen', opacity=0.8),
                transparentColor('tomato3', opacity=0.8),
                transparentColor('tomato', opacity=0.8))

    # k index for easy plotting
    rs  <-  unique(data$r)
    ks  <-  unique(data$k)

    # Set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

##  Row 1: 
    # Panel 1: r = 0.0
    dat  <-  subset(data, data$r == rs[1])
        par(omi=c(0.5, 0.5, 0.6, 0.5), mar = c(2.6,2.6,2.6,2.6), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(min(data$diffFemales),max(data$diffFemales)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different 
        abline(h=0, lwd=2, col='black')
        lines(diffFemales[k==ks[1]] ~ C[k==ks[1]], lwd=3, col=COLS[1], cex=1, data=dat)
        lines(diffFemales[k==ks[2]] ~ C[k==ks[2]], lwd=3, col=COLS[2], cex=1, data=dat)
        lines(diffFemales[k==ks[3]] ~ C[k==ks[3]], lwd=3, col=COLS[3], cex=1, data=dat)
        lines(diffFemales[k==ks[4]] ~ C[k==ks[4]], lwd=3, col=COLS[4], cex=1, data=dat)
        lines(diffFemales[k==ks[5]] ~ C[k==ks[5]], lwd=3, col=COLS[5], cex=1, data=dat)
        # axes
        axis(1, las=1,labels=NA)
        axis(2, las=1)
        proportionalLabel(1.2, 1.25, expression(paste("Gynodioecy")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(r)," = 0.0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35, 0.5, expression(paste(Delta," Female frequency")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
    rm(dat)

    # Panel 2: r = 0.05
    dat  <-  subset(data, data$r == rs[2])
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(min(data$diffFemales),max(data$diffFemales)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different values of k
        abline(h=0, lwd=2, col='black')
        lines(diffFemales[k==ks[1]] ~ C[k==ks[1]], lwd=3, col=COLS[1], cex=1, data=dat)
        lines(diffFemales[k==ks[2]] ~ C[k==ks[2]], lwd=3, col=COLS[2], cex=1, data=dat)
        lines(diffFemales[k==ks[3]] ~ C[k==ks[3]], lwd=3, col=COLS[3], cex=1, data=dat)
        # axes
        axis(1, las=1,labels=NA)
        axis(2, las=1,labels=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(r)," = ", 0.005)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # legend
        legend(
              x       =  usr[2]*0.975,
              y       =  usr[4],
              legend  =  c(
                          expression(paste(italic(k)~"="~italic(hat(k))%*%1.1)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.975)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.95)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.925)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.9))),
              lty     =  1,
              lwd     =  2,
              col     =  c(COLS[1],COLS[2],COLS[3],COLS[4],COLS[5]),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
    rm(dat)
    rm(rs)
    rm(ks)

##  Row 2: 
    rs  <-  unique(data2$r)
    ks  <-  unique(data2$k)
    # Panel 3: r = 0.0
    dat  <-  subset(data2, data2$r == rs[1])
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(min(data2$diffMales),max(data2$diffMales)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different values of k
        abline(h=0, lwd=2, col='black')
        lines(diffMales[k==ks[1]] ~ C[k==ks[1]], lwd=3, col=COLS[1], cex=1, data=dat)
        lines(diffMales[k==ks[2]] ~ C[k==ks[2]], lwd=3, col=COLS[2], cex=1, data=dat)
        lines(diffMales[k==ks[3]] ~ C[k==ks[3]], lwd=3, col=COLS[3], cex=1, data=dat)
        lines(diffMales[k==ks[4]] ~ C[k==ks[4]], lwd=3, col=COLS[4], cex=1, data=dat)
        lines(diffMales[k==ks[5]] ~ C[k==ks[5]], lwd=3, col=COLS[5], cex=1, data=dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(1.2, 1.25, expression(paste("Androdioecy")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(r)," = 0.0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35, 0.5, expression(paste(Delta," Male frequency")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.3, expression(paste(italic(C))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
    rm(dat)

    # Panel 4: r = 0.1
    dat  <-  subset(data2, data2$r == rs[2])
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(min(data2$diffMales),max(data2$diffMales)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different values of k
        abline(h=0, lwd=2, col='black')
        lines(diffMales[k==ks[1]] ~ C[k==ks[1]], lwd=3, col=COLS[1], cex=1, data=dat)
        lines(diffMales[k==ks[2]] ~ C[k==ks[2]], lwd=3, col=COLS[2], cex=1, data=dat)
        lines(diffMales[k==ks[3]] ~ C[k==ks[3]], lwd=3, col=COLS[3], cex=1, data=dat)
        lines(diffMales[k==ks[4]] ~ C[k==ks[4]], lwd=3, col=COLS[4], cex=1, data=dat)
        lines(diffMales[k==ks[5]] ~ C[k==ks[5]], lwd=3, col=COLS[5], cex=1, data=dat)
        # axes
        axis(1, las=1)
        axis(2, labels=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(r)," = ", 0.005)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(C))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
    rm(dat)
}


#' Multipanel plot summarizing results recursion simulations
#'
#' @title Figure 3 alternate
#' @author Colin Olito.
#' @export
Fig3Alt  <-  function() {

    # Import data
    data   <-  read.csv(file="./output/data/simResults/gyn-recess_dStar0.8_a1_sm0.1_add.csv", header=TRUE)
    data2  <-  read.csv(file="./output/data/simResults/and-recess_dStar0.8_a1_sm0.1_add.csv", header=TRUE)

    # Calculate equilibrium frequencies of M2, females, A, a
    data$q.m2        <-  (data$'F.12'/2) + (data$'F.14'/2) + data$'F.22' + (data$'F.23'/2) + data$'F.24' + (data$'F.34'/2) + data$'F.44' + 
                         (data$'G.12'/2) + (data$'G.14'/2) + data$'G.22' + (data$'G.23'/2) + data$'G.24' + (data$'G.34'/2) + data$'G.44'
    data$females     <-  data$'F.22' + data$'F.24' + data$'F.44' + data$'G.22' + data$'G.24' + data$'G.44'
    data$p.A         <-  data$'F.11' + data$'F.12' + (data$'F.13'/2) + (data$'F.14'/2) + data$'F.22' + (data$'F.23'/2) + (data$'F.24'/2) + 
                         data$'G.11' + data$'G.12' + (data$'G.13'/2) + (data$'G.14'/2) + data$'G.22' + (data$'G.23'/2) + (data$'G.24'/2)
    data$q.a         <-  (data$'F.13'/2) + (data$'F.14'/2) + (data$'F.23'/2) + (data$'F.24'/2) + data$'F.33' + data$'F.34' + data$'F.44' + 
                         (data$'G.13'/2) + (data$'G.14'/2) + (data$'G.23'/2) + (data$'G.24'/2) + data$'G.33' + data$'G.34' + data$'G.44' 
    data$diffFemales  <-  (data$females - data$ZHat)

   # Calculate equilibrium frequencies of M2, males, A, a
    data2$q.m2        <-  (data2$'F.12'/2) + (data2$'F.14'/2) + data2$'F.22' + (data2$'F.23'/2) + data2$'F.24' + (data2$'F.34'/2) + data2$'F.44' + 
                         (data2$'G.12'/2) + (data2$'G.14'/2) + data2$'G.22' + (data2$'G.23'/2) + data2$'G.24' + (data2$'G.34'/2) + data2$'G.44'
    data2$males       <-  data2$'F.22' + data2$'F.24' + data2$'F.44' + data2$'G.22' + data2$'G.24' + data2$'G.44'
    data2$p.A         <-  data2$'F.11' + data2$'F.12' + (data2$'F.13'/2) + (data2$'F.14'/2) + data2$'F.22' + (data2$'F.23'/2) + (data2$'F.24'/2) + 
                         data2$'G.11' + data2$'G.12' + (data2$'G.13'/2) + (data2$'G.14'/2) + data2$'G.22' + (data2$'G.23'/2) + (data2$'G.24'/2)
    data2$q.a         <-  (data2$'F.13'/2) + (data2$'F.14'/2) + (data2$'F.23'/2) + (data2$'F.24'/2) + data2$'F.33' + data2$'F.34' + data2$'F.44' + 
                         (data2$'G.13'/2) + (data2$'G.14'/2) + (data2$'G.23'/2) + (data2$'G.24'/2) + data2$'G.33' + data2$'G.34' + data2$'G.44' 
    data2$diffMales   <-  (data2$males - data2$ZHat)

    # Color scheme
    COLS  <-  c(transparentColor('dodgerblue', opacity=0.8),
              transparentColor('dodgerblue4', opacity=0.8),
              transparentColor('darkolivegreen', opacity=0.8),
              transparentColor('tomato3', opacity=0.8),
              transparentColor('tomato', opacity=0.8))

    # k index for easy plotting
    rs  <-  unique(data$r)
    ks  <-  unique(data$k)

    # Set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

##  Row 1: GYNODIOECY
    # Panel 1: r = 0.0
    dat  <-  subset(data, data$r == rs[1])
        par(omi=c(0.5, 0.5, 0.6, 0.5), mar = c(5,5,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(min(data$females),max(data$females)*1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different 
#        lines(ZHat[k==ks[1]] ~ C[k==ks[1]], lwd=2, col='black', cex=1, data=dat)
        lines(females[k==ks[1]] ~ C[k==ks[1]], lwd=2, col=COLS[1], cex=1, data=dat)
        lines(females[k==ks[2]] ~ C[k==ks[2]], lwd=2, col=COLS[2], cex=1, data=dat)
        lines(females[k==ks[3]] ~ C[k==ks[3]], lwd=2, col=COLS[3], cex=1, data=dat)
        lines(females[k==ks[4]] ~ C[k==ks[4]], lwd=2, col=COLS[4], cex=1, data=dat)
        lines(females[k==ks[5]] ~ C[k==ks[5]], lwd=2, col=COLS[5], cex=1, data=dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(1.2, 1.25, expression(paste("Gynodioecy")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.45, 1.1, expression(paste(italic(r)," =")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.635, 1.11, substitute(r,list(r=rounded(rs[1],precision=1))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35, 0.5, expression(paste("Eq. female frequency")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
    rm(dat)

    # Panel 2: r = 0.005
    dat  <-  subset(data, data$r == rs[2])
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(0,(max(data$females)*1.1)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different 
#        lines(ZHat[k==ks[1]] ~ C[k==ks[1]], lwd=2, col='black', cex=1, data=dat)
        lines(females[k==ks[1]] ~ C[k==ks[1]], lwd=2, col=COLS[1], cex=1, data=dat)
        lines(females[k==ks[2]] ~ C[k==ks[2]], lwd=2, col=COLS[2], cex=1, data=dat)
        lines(females[k==ks[3]] ~ C[k==ks[3]], lwd=2, col=COLS[3], cex=1, data=dat)
        lines(females[k==ks[4]] ~ C[k==ks[4]], lwd=2, col=COLS[4], cex=1, data=dat)
        lines(females[k==ks[5]] ~ C[k==ks[5]], lwd=2, col=COLS[5], cex=1, data=dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.45, 1.1, expression(paste(italic(r)," =")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.675, 1.11, substitute(r,list(r=rounded(rs[2],precision=3))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2],
              y       =  usr[4]*1.025,
              legend  =  c(
                          expression(paste(italic(k)~"="~italic(hat(k))%*%1.1)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.975)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.95)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.925)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.9))),
              lty     =  1,
              lwd     =  2,
              col     =  c(COLS[1],COLS[2],COLS[3],COLS[4],COLS[5]),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
    rm(dat)

##  Row 2: ANDRODIOECY
    rs  <-  unique(data2$r)
    ks  <-  unique(data2$k)

    # Panel 3: r = 0.0
    dat  <-  subset(data2, data2$r == rs[1])
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(0,max(data2$males)*1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different values of k
        lines(males[k==ks[1]] ~ C[k==ks[1]], lwd=2, col=COLS[1], cex=1, data=dat)
        lines(males[k==ks[2]] ~ C[k==ks[2]], lwd=2, col=COLS[2], cex=1, data=dat)
        lines(males[k==ks[3]] ~ C[k==ks[3]], lwd=2, col=COLS[3], cex=1, data=dat)
        lines(males[k==ks[4]] ~ C[k==ks[4]], lwd=2, col=COLS[4], cex=1, data=dat)
        lines(males[k==ks[5]] ~ C[k==ks[5]], lwd=2, col=COLS[5], cex=1, data=dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(1.2, 1.25, expression(paste("Androdioecy")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.1, expression(paste(italic(r)," = 0.0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35, 0.5, expression(paste("Eq. male frequency")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.3, expression(paste(italic(C))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
    rm(dat)

    # Panel 4: r = 0.1
    dat  <-  subset(data2, data2$r == rs[2])
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(0,max(data2$males)*1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different values of k
        lines(males[k==ks[1]] ~ C[k==ks[1]], lwd=2, col=COLS[1], cex=1, data=dat)
        lines(males[k==ks[2]] ~ C[k==ks[2]], lwd=2, col=COLS[2], cex=1, data=dat)
        lines(males[k==ks[3]] ~ C[k==ks[3]], lwd=2, col=COLS[3], cex=1, data=dat)
        lines(males[k==ks[4]] ~ C[k==ks[4]], lwd=2, col=COLS[4], cex=1, data=dat)
        lines(males[k==ks[5]] ~ C[k==ks[5]], lwd=2, col=COLS[5], cex=1, data=dat)
        # axes
        axis(1, las=1)
        axis(2, labels=NA)
#        proportionalLabel(0.5, 1.1, expression(paste(italic(r)," = ", 0.005)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(C))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
    rm(dat)
}


##############################################################
##############################################################
##  Exploratory plots

#' Exploratory plots: Invasion of dominant sterility allele into populations at 1-locus SA equilibrium. 
#'                    -- Additive fitness effects at SA locus
#' @title Invasion of dominant male sterility allele into populations
#' @author Colin Olito
#' @export
EQInv.Add  <-  function(df="./output/data/EQInvAnalyses/Gyn-partSelf-C25-delta80-strgSel-Add-EQInv.csv", wkSel=FALSE) {

    # Import data
    data  <-  read.csv(df, header=TRUE)

    # Color scheme
    COLS  <-  transparentColor('dodgerblue', opacity=0.2)

    # set plot limits for weak or strong selection
    if(wkSel)
        lims  <-  c(0,0.5)
    else(lims  <-  c(0,1))

    # Calculate 1-locus SA invasion criteria to illustrate 
    # boundaries for polymorphic populations
    if(any(colnames(data) == "C")) {
        sms  <-  seq(lims[1], lims[2], by=0.01)
        Ainv  <-  Inv.A.add(sms, C = data$C[1], delta = data$delta[1])
        Ainv[Ainv > 1]  <-  1.00001
        ainv  <-  Inv.a.add(sms, C = data$C[1], delta = data$delta[1])
    } else {
        sms  <-  seq(lims[1], lims[2], by=0.01)
        Ainv  <-  Inv.A.add(sms, C = 0, delta = 0)
        Ainv[Ainv > 1]  <-  1.00001
        ainv  <-  Inv.a.add(sms, C = 0, delta = 0)
    }

    # k index for easy plotting
    ks  <-  unique(data$k)
    rs  <-  unique(data$r)

    # Set plot layout
    layout.mat <- matrix(c(1:20), nrow=5, ncol=4, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

##  Row 1: r = 0
    ##  Panel 1: k = kCrit
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==rs[1] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==rs[1]]), precision=3)
        # Make plot
        par(omi=rep(0.4, 4), mar = c(3,3,0.75,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[1] & r==rs[1] & DiffEQInvEig != 0] ~ sm[k==ks[1] & r==rs[1] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.5, 1.25, expression(paste(hat(italic(k))%*%1.1," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.88, 1.23, substitute(k,list(k=rounded(ks[1],precision=2))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.6, 0.65, substitute(r,list(r=rs[1])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 2: k = kCrit - 0.1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==rs[1] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==rs[1]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==rs[1] & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==rs[1] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(hat(italic(k))%*%0.99," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.88, 1.23, substitute(k,list(k=rounded(ks[2],precision=2))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 3: k = kCrit - 0.2
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==rs[1] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==rs[1]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==rs[1] & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==rs[1] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(hat(italic(k))%*%0.95," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.88, 1.23, substitute(k,list(k=rounded(ks[3],precision=2))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 4: k = kCrit - 0.3
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[4] & data$r==rs[1] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[4] & data$r==rs[1]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[4] & r==rs[1] & DiffEQInvEig != 0] ~ sm[k==ks[4] & r==rs[1] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(hat(italic(k))%*%0.90," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.88, 1.23, substitute(k,list(k=rounded(ks[4],precision=2))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)


##  Row 2: r = 0.01
    ##  Panel 5: k = kCrit
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==rs[3] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==rs[3]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[1] & r==rs[3] & DiffEQInvEig != 0] ~ sm[k==ks[1] & r==rs[3] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.6, 0.7, substitute(r,list(r=rs[3])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, expression(paste(bold(E))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 6: k = kCrit - 0.1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==rs[3] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==rs[3]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==rs[3] & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==rs[3] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(F))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 7: k = kCrit - 0.2
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==rs[3] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==rs[3]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==rs[3] & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==rs[3] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(G))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 8: k = kCrit - 0.3
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[4] & data$r==rs[3] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[4] & data$r==rs[3]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[4] & r==rs[3] & DiffEQInvEig != 0] ~ sm[k==ks[4] & r==rs[3] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(H))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

##  Row 3: r = 0.02
    ##  Panel 9: k = kCrit
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==rs[4] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==rs[4]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[1] & r==rs[4] & DiffEQInvEig != 0] ~ sm[k==ks[1] & r==rs[4] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.6, 0.7, substitute(r,list(r=rs[4])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, expression(paste(bold(I))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 10: k = kCrit - 0.1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==rs[4] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==rs[4]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==rs[4] & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==rs[4] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(J))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 11: k = kCrit - 0.2
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==rs[4] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==rs[4]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==rs[4] & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==rs[4] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(K))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 12: k = kCrit - 0.3
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[4] & data$r==rs[4] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[4] & data$r==rs[4]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[4] & r==rs[4] & DiffEQInvEig != 0] ~ sm[k==ks[4] & r==rs[4] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(L))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)


##  Row 4: r = 0.1
    ##  Panel 13: k = kCrit
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==rs[5] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==rs[5]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[1] & r==rs[5] & DiffEQInvEig != 0] ~ sm[k==ks[1] & r==rs[5] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.6, 0.7, substitute(r,list(r=rs[5])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, expression(paste(bold(M))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 14: k = kCrit - 0.1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==rs[5] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==rs[5]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==rs[5] & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==rs[5] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(N))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 15: k = kCrit - 0.2
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==rs[5] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==rs[5]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==rs[5] & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==rs[5] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(O))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 16: k = kCrit - 0.3
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[4] & data$r==rs[5] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[4] & data$r==rs[5]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[4] & r==rs[5] & DiffEQInvEig != 0] ~ sm[k==ks[4] & r==rs[5] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(P))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

##  Row 5: r = 0.5
    ##  Panel 17: k = kCrit
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[1] & data$r==rs[6] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[1] & data$r==rs[6]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[1] & r==rs[6] & DiffEQInvEig != 0] ~ sm[k==ks[1] & r==rs[6] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.6, 0.7, substitute(r,list(r=rs[6])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(Q))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 18: k = kCrit - 0.1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[2] & data$r==rs[6] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[2] & data$r==rs[6]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[2] & r==rs[6] & DiffEQInvEig != 0] ~ sm[k==ks[2] & r==rs[6] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(R))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 19: k = kCrit - 0.2
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[3] & data$r==rs[6] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[3] & data$r==rs[6]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[3] & r==rs[6] & DiffEQInvEig != 0] ~ sm[k==ks[3] & r==rs[6] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(S))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

    ##  Panel 20: k = kCrit - 0.3
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==ks[4] & data$r==rs[6] & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==ks[4] & data$r==rs[6]]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = lims, ylim = lims, ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==ks[4] & r==rs[6] & DiffEQInvEig != 0] ~ sm[k==ks[4] & r==rs[6] & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
        # Overlay 1-locus invasion criteria
        lines(Ainv[Ainv<=1] ~ sms[Ainv<=1], lwd=2, col='black')
        lines(ainv ~ sms, lwd=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(T))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute("Pr(inv) ="~p, list(p = pInv)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        rm(pInv)

}


#' Exploratory plot: Illustration of different Selfing -- Inbreeding Depression relations. 
#' @title Comparing different C ~ delta functions
#' @author Colin Olito
#' @export
#' Multipanel plot comparing 2-locus recursion simulation results
#' with 1-locus predictions
#'
#' @title Figure 2 alternate
#' @author Colin Olito.
#' @export
compareCDelta  <-  function() {

# Calculate Different values of delta given deltaStar = 0.8 (as in the main text
C  <-  seq(0, 0.9, by=0.005)
linDelta     <-  deltaC(dStar = 0.8, C = C, a = 1, b = 0.5)
nonLinDelta  <-  deltaC(dStar = 0.8, C = C, a = 0.2, b = 0.5)

invGynLin     <-  invGyn(C = C, delta = linDelta) * 1.1
invGynNonLin  <-  invGyn(C = C, delta = nonLinDelta) * 1.1

Zhat.gyn.lin  <-  c()
Zhat.gyn.nonlin  <-  c()
    for(i in 1:length(C)) {
        Zhat.gyn.lin[i]     <-  Zhat.gyn(k = invGynLin[i], C = C[i], delta = linDelta[i])
        Zhat.gyn.nonlin[i]  <-  Zhat.gyn(k = invGynNonLin[i], C = C[i], delta = nonLinDelta[i])
    } 


    # Set plot layout
    layout.mat <- matrix(c(1:2), nrow=1, ncol=2, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

### Equilibrium Female Frequencies
##  Panel 1: 
        par(omi=rep(0.3, 4), mar = c(3.5,3.5,2,2), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different 
        lines(linDelta ~ C, lwd=2, col='black', cex=1)
        lines(nonLinDelta ~ C, lwd=2, lty=2, col='black', cex=1)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.05, 1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(C))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35, 0.5, expression(paste(delta)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        #legend
        legend(
              x       =  usr[2]*0.99,
              y       =  usr[4],
              legend  =  c(
                          expression(paste("Linear: ",italic(a)," = ",1)),
                          expression(paste("Non-linear: ",italic(a)," = ",0.2))),
              lty     =  c(1,2),
              lwd     =  2,
              col     =  'black',
              cex     =  0.75,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )


    # Import data
    data   <-  read.csv(file="./output/data/simResults/gyn-recess_dStar0.8_a1_sm0.4_add.csv", header=TRUE)
    data2  <-  read.csv(file="./output/data/simResults/gyn-recess_dStar0.8_a0.2_sm0.4_add.csv", header=TRUE)

    # Calculate equilibrium frequencies of M2, females, A, a
    data$q.m2         <-  (data$'F.12'/2) + (data$'F.14'/2) + data$'F.22' + (data$'F.23'/2) + data$'F.24' + (data$'F.34'/2) + data$'F.44' + 
                          (data$'G.12'/2) + (data$'G.14'/2) + data$'G.22' + (data$'G.23'/2) + data$'G.24' + (data$'G.34'/2) + data$'G.44'
    data$females      <-  data$'F.22' + data$'F.24' + data$'F.44' + data$'G.22' + data$'G.24' + data$'G.44'
    data$p.A          <-  data$'F.11' + data$'F.12' + (data$'F.13'/2) + (data$'F.14'/2) + data$'F.22' + (data$'F.23'/2) + (data$'F.24'/2) + 
                          data$'G.11' + data$'G.12' + (data$'G.13'/2) + (data$'G.14'/2) + data$'G.22' + (data$'G.23'/2) + (data$'G.24'/2)
    data$q.a          <-  (data$'F.13'/2) + (data$'F.14'/2) + (data$'F.23'/2) + (data$'F.24'/2) + data$'F.33' + data$'F.34' + data$'F.44' + 
                          (data$'G.13'/2) + (data$'G.14'/2) + (data$'G.23'/2) + (data$'G.24'/2) + data$'G.33' + data$'G.34' + data$'G.44' 
    data$diffFemales  <-  (data$females - data$ZHat)

   # Calculate equilibrium frequencies of M2, males, A, a
    data2$q.m2         <-  (data2$'F.12'/2) + (data2$'F.14'/2) + data2$'F.22' + (data2$'F.23'/2) + data2$'F.24' + (data2$'F.34'/2) + data2$'F.44' + 
                           (data2$'G.12'/2) + (data2$'G.14'/2) + data2$'G.22' + (data2$'G.23'/2) + data2$'G.24' + (data2$'G.34'/2) + data2$'G.44'
    data2$females      <-  data2$'F.22' + data2$'F.24' + data2$'F.44' + data2$'G.22' + data2$'G.24' + data2$'G.44'
    data2$p.A          <-  data2$'F.11' + data2$'F.12' + (data2$'F.13'/2) + (data2$'F.14'/2) + data2$'F.22' + (data2$'F.23'/2) + (data2$'F.24'/2) + 
                           data2$'G.11' + data2$'G.12' + (data2$'G.13'/2) + (data2$'G.14'/2) + data2$'G.22' + (data2$'G.23'/2) + (data2$'G.24'/2)
    data2$q.a          <-  (data2$'F.13'/2) + (data2$'F.14'/2) + (data2$'F.23'/2) + (data2$'F.24'/2) + data2$'F.33' + data2$'F.34' + data2$'F.44' + 
                           (data2$'G.13'/2) + (data2$'G.14'/2) + (data2$'G.23'/2) + (data2$'G.24'/2) + data2$'G.33' + data2$'G.34' + data2$'G.44' 
    data2$diffFemales  <-  (data2$females - data2$ZHat)

    # Color scheme
    COLS  <-  c(transparentColor('dodgerblue', opacity=0.8),
                transparentColor('dodgerblue4', opacity=0.8),
                transparentColor('darkolivegreen', opacity=0.8),
                transparentColor('tomato3', opacity=0.8),
                transparentColor('tomato', opacity=0.8))

    # k index for easy plotting
    rs  <-  unique(data$r)
    ks  <-  unique(data2$k)

    dat1  <-  subset(data, data$k == ks[1])
    dat2  <-  subset(data2, data2$k == ks[1])

        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(0,0.22), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different 
        lines(Zhat.gyn.lin ~ C, lwd=2, col='black')
        lines(Zhat.gyn.nonlin ~ C, lwd=2, lty=2, col='black')
        lines(females[r==rs[1]] ~ C[r==rs[1]], lwd=2, col='grey70', lty=1, cex=1, data=dat1)
        lines(females[r==rs[1]] ~ C[r==rs[1]], lwd=2, col='grey70', lty=2, cex=1, data=dat2)
        lines(females[r==rs[2]] ~ C[r==rs[2]], lwd=2, col='grey60', lty=1, cex=1, data=dat1)
        lines(females[r==rs[2]] ~ C[r==rs[2]], lwd=2, col='grey60', lty=2, cex=1, data=dat2)
        lines(females[r==rs[3]] ~ C[r==rs[3]], lwd=2, col='grey50', lty=1, cex=1, data=dat1)
        lines(females[r==rs[3]] ~ C[r==rs[3]], lwd=2, col='grey50', lty=2, cex=1, data=dat2)
        lines(females[r==rs[4]] ~ C[r==rs[4]], lwd=2, col='grey40', lty=1, cex=1, data=dat1)
        lines(females[r==rs[4]] ~ C[r==rs[4]], lwd=2, col='grey40', lty=2, cex=1, data=dat2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.35, 0.5, expression(paste("Eq. female frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.3, expression(paste(italic(C))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*0.99,
              y       =  usr[4],
              legend  =  c(
                          expression(paste("Linear: ",italic(a)," = ",1)),
                          expression(paste("Non-linear: ",italic(a)," = ",0.2))),
              lty     =  c(1,2),
              lwd     =  2,
              col     =  'black',
              cex     =  0.75,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
    proportionalLabel(0.115, 0.93, expression(paste(italic(r)~"="~0.0)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.165, 0.68, expression(paste(italic(r)~"="~0.005)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.15, 0.58, expression(paste(italic(r)~"="~0.01)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.15, 0.42, expression(paste(italic(r)~"="~0.05)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.15, 0.26, expression(paste("1 locus")), cex=0.75, adj=c(0.5, 0.5), xpd=NA)

}