#--------------------------------------------------
# Data preparation.  For 2-state mcmc and maps.
# Keep good states and tips.  Compute sampling.f.
#--------------------------------------------------

state.names <- c(focal.state, "Other")
sn <- c(t4[focal.state], Other = "O")

dat <- read.csv(paste("../../../data/", genus, ".csv", sep=""), comment.char="#")

# Arrange data structure to compute sampling.f.  Remove empty states.
dat <- data.frame(Name=dat$Name, State=dat$SexSyst4, OnTree=dat$OnTree)
i.drop <- which(dat$State == "")
if (length(i.drop) > 0)
    dat <- dat[-i.drop,]
# not focal.state at all
dat$State[!unlist(lapply(strsplit(dat$State, split="|", fixed=T),
                           function(x) focal.state %in% x))] <- "Other"
dat$State[!(dat$State %in% c(focal.state, "Other"))] <- 
                                    paste(focal.state, "Other", sep="|")

samp.f <- get.sampf(dat)
samp.f <- samp.f[state.names]

# Arrange states for bisse.  Must be on tree.  Can be uncertain.
dat <- subset(dat, OnTree == "yes")
states <- structure(rep(NA, nrow(dat)), names=dat$Name)
states[dat$State == state.names[1]] <- 0
states[dat$State == state.names[2]] <- 1

# Check if it's worth proceeding
if (any(table(states) < 2) | length(table(states)) < 2)
    stop("Not enough species in each state to proceed.")

# Phylogenies
phy.all <- read.tree(paste("../../../data/", genus, ".tre", sep=""))
dropme <- phy.all[[1]]$tip.label[!(phy.all[[1]]$tip.label %in% names(states))]

# Misc
outdir <- "."
if (!file.exists(outdir))
    dir.create(outdir)

if (!exists("i.phy"))
    i.phy <- seq_along(phy.all)
