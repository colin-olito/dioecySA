###############
# DEPENDENCIES
###############
# library(...)


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






##############################################################
##############################################################
##  Exploratory figures



#' Exploratory plots: Invasion of dominant male sterility allele into populations
#'                    at 1-locus SA equilibrium. 
#'                    -- Obligate outcrossing
#'                    -- Additive fitness effects at SA locus
#' @title Fig.E1: Invasion of dominant male sterility allele into populations
#'                    at 1-locus SA equilibrium. 
#'                    -- Obligate outcrossing
#'                    -- Additive fitness effects at SA locus
#' @export
EQInv.ObOut.Add.Gyno.Dom  <-  function() {

    # Import data
    data  <-  read.csv("./output/data/EQInvAnalyses/ObOut-Add-EQInv.csv", header=TRUE)

    # Color scheme
    COLS  <-  transparentColor('dodgerblue4', opacity=0.2)

    # Calculate 1-locus SA invasion criteria to illustrate 
    # boundaries for polymorphic populations
    sms  <-  seq(0,1,by=0.01)
    Ainv  <-  inv.AKid(sms)
    Ainv[Ainv > 1]  <-  1.00001
    ainv  <-  inv.aKid(sms)

    # Set plot layout
    layout.mat <- matrix(c(1:15), nrow=5, ncol=3, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

##  Row 1: r = 0
    ##  Panel 1: k = 1
        # Calculate proportion parameter space where sterility allele can invade
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==1 & data$r==0 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==1 & data$r==0]), precision=3)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==0.875 & data$r==0 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==0.875 & data$r==0]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==0.875 & r==0 & DiffEQInvEig != 0] ~ sm[k==0.875 & r==0 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==0.75 & data$r==0 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==0.75 & data$r==0]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==0.75 & r==0 & DiffEQInvEig != 0] ~ sm[k==0.75 & r==0 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==1 & data$r==0.01 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==1 & data$r==0.01]), precision=3)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==0.875 & data$r==0.01 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==0.875 & data$r==0.01]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==0.875 & r==0.01 & DiffEQInvEig != 0] ~ sm[k==0.875 & r==0.01 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==0.75 & data$r==0.01 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==0.75 & data$r==0.01]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==0.75 & r==0.01 & DiffEQInvEig != 0] ~ sm[k==0.75 & r==0.01 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==1 & data$r==0.02 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==1 & data$r==0.02]), precision=3)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==0.875 & data$r==0.02 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==0.875 & data$r==0.02]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==0.875 & r==0.02 & DiffEQInvEig != 0] ~ sm[k==0.875 & r==0.02 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==0.75 & data$r==0.02 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==0.75 & data$r==0.02]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==0.75 & r==0.02 & DiffEQInvEig != 0] ~ sm[k==0.75 & r==0.02 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==1 & data$r==0.1 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==1 & data$r==0.1]), precision=3)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==0.875 & data$r==0.1 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==0.875 & data$r==0.1]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==0.875 & r==0.1 & DiffEQInvEig != 0] ~ sm[k==0.875 & r==0.1 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==0.75 & data$r==0.1 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==0.75 & data$r==0.1]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==0.75 & r==0.1 & DiffEQInvEig != 0] ~ sm[k==0.75 & r==0.1 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==1 & data$r==0.5 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==1 & data$r==0.5]), precision=3)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==0.875 & data$r==0.5 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==0.875 & data$r==0.5]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==0.875 & r==0.5 & DiffEQInvEig != 0] ~ sm[k==0.875 & r==0.5 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
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
            pInv  <-  rounded(length(data$DiffEQInvEig[data$k==0.75 & data$r==0.5 & data$DiffEQInvEig != 0]) / 
                              length(data$DiffEQInvEig[data$k==0.75 & data$r==0.5]), precision=3)
        # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(sf[k==0.75 & r==0.5 & DiffEQInvEig != 0] ~ sm[k==0.75 & r==0.5 & DiffEQInvEig != 0], pch=21, col=NA, cex=1, bg=COLS, data=data)
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