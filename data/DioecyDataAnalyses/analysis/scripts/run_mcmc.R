#--------------------------------------------------
# Fit BiSSE or MuSSE model, constrained to state-independent diversification.
# Correct for biased sampling.  Include polymorphic tips.
#
# Call this script from within the directory of a genus.
#--------------------------------------------------

args <- commandArgs(TRUE)
if (length(args) != 1)
    stop("Need to provide focal state [Herm, Dio, xDio, xMono; or 4]")
focal.state <- args[1]

genus <- rev(strsplit(getwd(), "/")[[1]])[1]

source("../../scripts/commonP.R")
source("../../scripts/common4.R")
source("../../scripts/common2.R")

# Optional.  If commented-out, will run all trees.
# i.phy <- 31:34

# Empirical prior
prior.rates <- c("Herm" = 0.6, "Dio" = 0.4, "4" = 0.35)
prior.rate <- prior.rates[focal.state]

# Read in cleaned data
if (focal.state == "4")
{
    source("../../scripts/prep4.R")
} else {
    source("../../scripts/prep2.R")
}

#--------------------------------------------------
# Run MCMC on each tree
#--------------------------------------------------

run.mcmc <- function(i.tree)
{
    message(paste("Running model for", genus, "tree", i.tree))

    # Prepare the one tree.
    phy <- phy.all[[i.tree]]
    phy$edge.length <- phy$edge.length / max(branching.times(phy))
    age <- max(branching.times(phy))  # in case rescaling is abandoned
    phy <- drop.species(phy, dropme, quiet=T)

    # Prepare the Mk model
    if (focal.state == "4")
    {
        lnL.full <- make.musse(phy, states, k, sampling.f=samp.f)
        lnL <- constrain.to.mkn(lnL.full, state.names)
    } else {
        lnL.full <- make.bisse(phy, states, sampling.f=samp.f)
        lnL <- constrain.to.mkn(lnL.full, sn)
    }
    lnL.names <- argnames(lnL)
    q.names <- lnL.names[grep("q", lnL.names)]
    n.par <- length(lnL.names)

    # Prior
    lp <- 1/starting.point.bd(phy)[1] # prior rate for lambda
    prior <- make.prior.exponential(c(lp, lp/10, rep(1/prior.rate, length(q.names))))

    # Pilot: determine chain parameters, remove most burnin
    nsteps <- 250
    burnin <- 50
    upper <- age*100
    w <- age
    par.start <- rep(age/2, length(argnames(lnL)))
    ans <- mcmc(lnL, par.start, nsteps=nsteps, w=w, upper=upper,
                prior=prior, print.every=0)
    init <- suggest.mcmc.params(ans[-seq(burnin),])

    # Real chain
    nsteps <- 1000
    outfile <- paste(outdir, "/mcmc_tree", zfill(i.tree, 3), ".csv", sep="")
    mcmc(lnL, init$par.start, nsteps=nsteps, w=init$w,
         upper=init$upper, prior=prior, print.every=0, save.every=100,
         save.file=outfile)

    return(0)
}

junk <- mclapply(i.phy, run.mcmc, mc.cores=detectCores())
