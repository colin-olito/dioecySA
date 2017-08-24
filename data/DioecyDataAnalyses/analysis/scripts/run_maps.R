#--------------------------------------------------
# Run stochastic mapping for 2-state or 4-state trait.
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

if (focal.state != "4")
    source("../../scripts/common2.R")

# Optional.  If commented-out, will run all trees.
# i.phy <- 31:34

# Read in cleaned data
if (focal.state == "4")
{
    source("../../scripts/prep4.R")
    i.uncertain <- grep.full("|", states1, fixed=T)
} else {
    source("../../scripts/prep2.R")
    i.uncertain <- which(is.na(states))
}
# If any tips are of uncertain state, "states" will be modified later.

#--------------------------------------------------
# Run stochastic map
#--------------------------------------------------

message(paste("Starting stochastic mapping for", genus))

run.stoch <- function(i.tree)
{
    outfile <- paste(outdir, "/maps_tree", zfill(i.tree, 3), ".csv", sep="")

    # Read in the MCMC rate samples for that tree.
    infile <- paste(outdir, "/mcmc_tree", zfill(i.tree, 3), ".csv", sep="")
    if (!file.exists(infile))
    {
        message(paste("No file found:", infile))
        return(0)
    }
    dat <- read.csv(infile)
    q.names <- colnames(dat)[grep("q", colnames(dat))]

    # Prepare the tree.
    phy <- phy.all[[i.tree]]
    phy$edge.length <- phy$edge.length / max(branching.times(phy))
    phy <- drop.species(phy, dropme, quiet=T)

    # Assemble the matrix of transition rates to use for the mapping.
    # Could be simply all samples, or draw randomly (with replacement).
    # Will obtain one mapping per row.
    i.samp <- seq(nrow(dat))
    pars <- dat[i.samp, q.names]
    pars <- as.list(data.frame(t(pars))) # to allow mclapply()

    # Generate the stochastic character mappings.
    if (focal.state != "4")
    {
        if (length(i.uncertain) > 0)
            states[i.uncertain] <- sample(0:1, size=length(i.uncertain), replace=T)
        lik <- make.mk2(phy, states, control=list(method="mk2"))
        maps <- mclapply(pars, function(x) safe.asr.stoch(lik, x), 
                         mc.cores=detectCores())
    } else {
        # (If there are no uncertain tip states, can use a faster method.)
        if (length(i.uncertain) == 0)
        {
            lik <- make.mkn(phy, states, k=k)
            maps <- mclapply(pars, function(x) safe.asr.stoch(lik, x), 
                             mc.cores=detectCores())
        } else {
            maps <- mclapply(pars, function(x) safe.asr.stoch.uncertain(phy, 
                                               states1, state.names,
                                               i.uncertain, x),
                             mc.cores=detectCores())
        }
    }

    # Weed out any mapping attempts that failed.
    i.good <- which(!sapply(maps, is.null))
    maps <- maps[i.good]
    i.samp <- i.samp[i.good]

    # Summarize the results
    ans.trans <- t(sapply(maps, summarize.transitions, sn, flatten=T))
    ans.root <- sapply(maps, get.root.state, sn)
    ans.time <- t(sapply(maps, summarize.durations, phy, sn))

    # Record the results
    ans2 <- data.frame("i" = i.samp, ans.trans, root=ans.root, ans.time)
    write.csv0(ans2, filename=outfile)

    return(0)
}

junk <- mclapply(i.phy, run.stoch, mc.cores=detectCores())
