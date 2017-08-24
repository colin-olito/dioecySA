##############################################################
### Histograms of dXY and nXY - nYX, colored by root state ###
##############################################################

library("ggplot2")
library("gridExtra")
library("grid")

col.state <- c(H = "gold1", M = "green4", G = "blue2", D = "red4", O = "gray")

#--------------------------------------------------
# 1. Read in and bin results
#    For all transitions
#--------------------------------------------------

# See Fig4-S10-14.R to obtain:
#   dat.prop.real
#   XY.YX.all
# But use this instead of the one in that file
get.prop <- function(dat, XY.YX)
{
    # extract the pathway of interest
    dat.leg <- mclapply(dat, get.maps.pathway, XY.YX)
    dat.leg <- dat.leg[!is.na(dat.leg)]
    dat <- do.call(rbind, dat.leg)

    # computations in maps
    dat$sum <- dat[[XY.YX[1]]] + dat[[XY.YX[2]]]
    dat <- dat[dat$sum > 0,]
    dat$prop <- dat[[XY.YX[1]]] / dat$sum 
    dat$diff <- dat[[XY.YX[1]]] - dat[[XY.YX[2]]] # new

    return(dat)
}

dat.all <- dat.prop.real
rm(dat.prop.real)

#--------------------------------------------------
# Plot a pathway
#--------------------------------------------------

options(stringsAsFactors = TRUE)

XY.all <- sapply(XY.YX.all, function(x) x[1])

xlim.num <- c(-15.5, 15.5)  # but does truncate for HM and HO

theme.bare <- theme(line = element_blank(),
                    axis.text.y = element_blank(),
                    title = element_blank(),
                    legend.position = "none",
                    panel.background = element_rect(fill = "white"),
                    plot.margin = rep(unit(0, "null"), 4)
                   )

plots.all <- list()
for (XY in XY.all)
{
    dat <- dat.all[[XY]]
    dat$root <- factor(dat$root, levels=rev(c("H", "G", "M", "D", "O")))

    ### dXY ###

    p.prop <- ggplot(dat, aes(x = prop)) +
          xlab("proportional asymmetry of transitions") +
          ylab("number of stochastic mappings") +
          ylim(0, 0.5) +
          theme(legend.text = element_text(face = "italic")) +
          scale_x_continuous(breaks=seq(0, 1, by=0.5)) +
          theme.bare

    p.prop.hist <- p.prop +
          geom_histogram(aes(y = (..count..)/sum(..count..)), 
                         position='stack', breaks=seq(0, 1, length.out=30))

    p.prop.root <- p.prop.hist + 
          aes(group = root, fill = root) + 
          scale_fill_manual(values=col.state)

    ### nXY - nYX ###

    p.num <- ggplot(dat, aes(x = diff)) +
          xlab("difference in number of transitions") +
          ylab("number of stochastic mappings") +
          ylim(0, 0.4) +
          theme(legend.text = element_text(face = "italic")) +
          geom_vline(xintercept = 0, color = "black", lty=3) +
          theme.bare

    p.num.hist <- p.num +
          geom_histogram(aes(y = (..count..)/sum(..count..)), position='stack', 
                   breaks = seq(xlim.num[1], xlim.num[2], by=1))

    p.num.root <- p.num.hist + 
          aes(group = root, fill = root) + 
          scale_fill_manual(values=col.state)

    plots.all[[XY]] <- arrangeGrob(p.prop.root, p.num.root, ncol=2, top=XY)
}

pdf(file="FigS4.pdf", width=5, height=10)
do.call(grid.arrange, c(plots.all, list(ncol=1)))
dev.off()
