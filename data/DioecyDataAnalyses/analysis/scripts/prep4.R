#--------------------------------------------------
# Data preparation.  For 4-state mcmc and maps.
# Keep good states and tips.  Compute sampling.f.
#--------------------------------------------------

# All species
dat <- read.csv(paste("../../../data/", genus, ".csv", sep=""), comment.char="#")
dat$State <- dat$SexSyst4

# Only species on tree
dat.phy <- subset(dat, OnTree == "yes")
ss <- structure(dat.phy$State, names=dat.phy$Name)

# Convert to matrix state structure
states <- as.multistate.matrix(ss, empty.action="drop", multi.action="equal")
state.names <- colnames(states)

# Actually, need non-matrix structure for uncertain tips
states1 <- ss
i <- which(states1 == "")
if (length(i) > 0)
    states1 <- states1[-i]

# Sampling
samp.f <- get.sampf(dat)
samp.f <- samp.f[state.names]

# Phylogenies
phy.all <- read.tree(paste("../../../data/", genus, ".tre", sep=""))
dropme <- phy.all[[1]]$tip.label[!(phy.all[[1]]$tip.label %in% rownames(states))] # if using drop in multistate above

# Misc
k <- length(state.names)  # number of states
sn <- t4[state.names]     # state name abbreviations

outdir <- "."
if (!file.exists(outdir))
    dir.create(outdir)

if (!exists("i.phy"))
    i.phy <- seq_along(phy.all)
