t4 <- c("H", "M", "G", "D")
names(t4) <- traits4

#--------------------------------------------------
# Utilities for model fitting
#--------------------------------------------------

### Convert vector of states to matrix of states (suitable for MuSSE).
# dat = character vector of states, with names = species names
#       multiple values separated by sep
# ans = species names as rownames, one column for each trait value
#       shows presence/absence [1/0] for that value
# empty.action = what to do with taxa lacking trait data
#     "none" = leave them alone (leave as 0)
#     "drop" = remove from trait matrix
#     "freq" = assign them the observed state frequencies
#     a number or vector = set all values to this
# multi.action = what to do with taxa in multiple states
#     "all" = code each state as 1
#     "equal" = code each of the n states as 1/n
#     "drop" = code each state as 0
# state.names = can be provided to ensure consistent columns even if a state
#               is absent
as.multistate.matrix <- function(dat, empty.action = "none", 
                                 multi.action = "all", sep = "|", 
                                 state.names = NULL)
{
    # Get info from the trait data.
    if (is.null(state.names))
        state.names <- sort(unique(unlist(strsplit(dat, sep, fixed=T))))
    species.names <- names(dat)

    # Prepare the matrix structure.
    ans <- as.data.frame(matrix(0, nrow=length(dat), ncol=length(state.names)))
    rownames(ans) <- species.names 
    names(ans) <- state.names

    # Fill in the "presence" values.
    for (sp in species.names)
    {
        # All values for this species
        st.all <- strsplit(dat[sp], sep, fixed=T)[[1]]

        # Skip empties
        i <- which(st.all %in% c("", NA))
        if (length(i) > 0)
            st.all <- st.all[-i]

        if (multi.action == "equal")
        {
            fill <- 1 / length(st.all)
        } else if (multi.action == "drop" & length(st.all) > 1) {
            fill <- 0
        } else {
            fill <- 1
        }

        for (st in st.all)
            ans[sp, st] <- fill
    }

    freq <- colSums(ans) / sum(ans)

    # Handle species with "absence" for all values.
    if (empty.action[1] != "none")
    {
        i <- which(rowSums(ans) == 0)
        if (length(i) > 0)
        {
            if (empty.action[1] == "drop")
            {
                ans <- ans[-i,]
            } else if (empty.action[1] == "freq")
            {
                for (i1 in i)
                    ans[i1,] <- freq
            } else {
                for (i1 in i)
                    ans[i1,] <- empty.action
            }
        }
    }

    return(as.matrix(ans))
}

# could reduce redundancy among the next four functions...

# Constrain MuSSE to Mkn
# sn = state names
constrain.to.mkn <- function(lik, sn)
{
    from0 <- "bisse" %in% class(lik)  # to adjust indices below
    an <- argnames(lik)

    # Remove state-dependent speciation and extinction
    lambda.names <- an[grep("lambda", an)]
    mu.names <- an[grep("mu", an)]

    # Rename the transition rates, to make comparisons across genera easier
    # (t4 is defined above)
    q.old <- an[grep("q", an)]
    q.new <- q.old
    for (i in seq_along(sn))
    {
        if (!from0)
        {
            q.new <- gsub(as.character(i), t4[sn[i]], q.new)
        } else {
            q.new <- gsub(as.character(i-1), sn[i], q.new)
        }
    }

    constr <- as.list(c(paste(lambda.names, "~ lambda"), 
                        paste(mu.names, "~ mu"),
                        paste(q.old, q.new, sep=" ~ ")
                        ))
    lik.ans <- constrain(lik, formulae=lapply(constr, formula),
                         extra=c("lambda", "mu", q.new))
    return(lik.ans)
}

# Constrain MuSSE to equal-transition-rates Mkn
constrain.to.mk1 <- function(lik)
{
    an <- argnames(lik)
    # Remove state-dependent speciation and extinction
    lambda.names <- an[grep("lambda", an)]
    mu.names <- an[grep("mu", an)]

    # Rename the transition rates, to make comparisons across genera easier
    # (t4 is defined above)
    q.old <- an[grep("q", an)]
    q.new <- "q"

    constr <- as.list(c(paste(lambda.names, "~ lambda"), 
                        paste(mu.names, "~ mu"),
                        paste(q.old, q.new, sep=" ~ ")
                        ))
    lik.ans <- constrain(lik, formulae=lapply(constr, formula),
                         extra=c("lambda", "mu", q.new))
    return(lik.ans)
}

rename.rates.sdd <- function(lik, sn)
{
    # Rename the diversification and transition rates, to make comparisons
    # across genera easier.  Only for BiSSE (not MuSSE).
    an <- argnames(lik)
    q.old <- an[grep("q", an)]
    lam.old <- an[grep("lambda", an)]
    mu.old <- an[grep("mu", an)]
    q.new <- c(paste("q", sn[1], sn[2], sep=""), 
               paste("q", sn[2], sn[1], sep=""))
    lam.new <- c(paste("lambda", sn[1], sep=""), 
                 paste("lambda", sn[2], sep=""))
    mu.new <- c(paste("mu", sn[1], sep=""), paste("mu", sn[2], sep=""))

    constr <- as.list(c(paste(lam.old, lam.new, sep=" ~ "),
                        paste(mu.old, mu.new, sep=" ~ "),
                        paste(q.old, q.new, sep=" ~ ")
                        ))
    lik.ans <- constrain(lik, formulae=lapply(constr, formula),
                         extra=c(lam.new, mu.new, q.new))
    return(lik.ans)
}

rename.rates <- function(lik, sn)
{
    # Rename the transition rates, to make comparisons across genera easier
    # (t4 is defined above)
    an <- argnames(lik)
    q.old <- an[grep("q", an)]
    q.new <- q.old
    for (i in seq_along(sn))
        q.new <- gsub(as.character(i), t4[sn[i]], q.new)

    constr <- as.list(paste(q.old, q.new, sep=" ~ "))
    lik.ans <- constrain(lik, formulae=lapply(constr, formula),
                         extra=q.new)
    return(lik.ans)
}

# All transition rates equal (and named "q")
rename.rates.equal <- function(lik)
{
    # Rename the transition rates, to make comparisons across genera easier
    an <- argnames(lik)
    q.old <- an[grep("q", an)]
    q.new <- "q"

    constr <- as.list(paste(q.old, q.new, sep=" ~ "))
    lik.ans <- constrain(lik, formulae=lapply(constr, formula),
                         extra=q.new)
    return(lik.ans)
}

# Compute sampling probabilities by state
# Requires columns named: Name, State, OnTree
# keep.na = if keeping tips that lack state data
get.sampf <- function(dat, keep.na = TRUE)
{
    # all
    st2 <- structure(dat$State, names=dat$Name)
    st2 <- as.multistate.matrix(st2, empty.action="none", multi.action="equal")

    # on the tree
    dat.phy <- subset(dat, OnTree == "yes")
    st1 <- structure(dat.phy$State, names=dat.phy$Name)
    st1 <- as.multistate.matrix(st1, empty.action="none", multi.action="equal")

    # states in genus but not on tree
    i <- which(!(colnames(st2) %in% colnames(st1)))
    if (length(i) > 0)
    {
        st1 <- cbind(st1, temp = 0)
        if (length(i) > 1) # could be more elegant...
            st1 <- cbind(st1, temp = 0)
        colnames(st1)[colnames(st1)=="temp"] <- colnames(st2)[i]
        st1 <- st1[,colnames(st2)]
    }
    stopifnot(all(colnames(st2) == colnames(st1)))

    N <- nrow(st2)                                 # total number of species
    x <- c(colSums(st1), xxx=sum(rowSums(st1)==0)) # on tree for each state
    z <- colSums(st2)                              # all known
    p <- z / sum(z)                                # prob of each state

    # Observed vs expected proportion of tips in each state on the tree
    n <- length(x)
    if (keep.na)
    {
        # unknowns are weighted by observed proportions
        samp.f <- (x[-n] + x[n]*p) / (N * p)
    } else {
        samp.f <- (x[-n]) / (sum(x[-n]) * p)
    }
    # If allowing keep.na = F, need to deal with samp.f > 1.

    return(samp.f)
}

# Given some mcmc samples, suggest control parameters
suggest.mcmc.params <- function(ans)
{
    i.par <- seq(from=2, to=ncol(ans)-1)

    w <- round(apply(ans[, i.par, drop=F], 2, max) - 
               apply(ans[,i.par, drop=F], 2, min), 2)
    upper <- round(apply(ans[, i.par, drop=F], 2, max) * 5, 2)
    par.start <- as.numeric(ans[which.max(ans$p), i.par])

    return(list(w=w, upper=upper, par.start=par.start))
}

#--------------------------------------------------
# Utilities for stochastic mapping
#--------------------------------------------------

safe.asr.stoch <- function(lik, pars, tries=100)
{
    for (i in seq(tries))
    {
        ans <- try(asr.stoch(lik, pars, n=1), silent=T)
        if (class(ans) == "history")
            return(ans)
    }
    # message("Failed to obtain a stochastic map")
    return(NULL)
}

# states1 has uncertain/multiple tip states for species numbers i.uncertain
safe.asr.stoch.uncertain <- function(phy, states1, state.names, i.uncertain, pars)
{
    k <- length(state.names)
    for (i in i.uncertain)
        states1[i] <- sample(strsplit(states1[i], split="|", fixed=T)[[1]], 1)
    statesN <- as.multistate.matrix(states1, state.names=state.names)
    lik <- make.mkn(phy, statesN, k=k)
    safe.asr.stoch(lik, pars)
}

## Given a history object (h, from simulation or stochastic mapping), 
## report how many transitions have happened.
## Not speed optimized!  But tested.
## sn = optional names for states (to over-ride the numbers in h)
## flatten = return a vector rather than a matrix
summarize.transitions <- function(h, sn = NULL, flatten = FALSE) {
  if (class(h) != "history")
    stop("Need to provide a history object to summarize.transitions()")

  # translation between state names and indices
  # (needed because states may start from 0 or 1
  nst <- length(h$states)
  sidx <- seq(nst)
  names(sidx) <- h$states

  # Use a matrix to keep track of how many of each kind of change.
  # (Diagonal elements will always be 0.)
  emptycount <- matrix(0, nrow=nst, ncol=nst)

  # Construct a count matrix for one edge
  count.edge <- function(edge) {
    count <- emptycount
    if (nrow(edge) > 1) {
      for (i in seq(nrow(edge)-1)) {
        # translate state names into indices of the "count" matrix
        i1 <- sidx[as.character(edge[i,2])]
        i2 <- sidx[as.character(edge[i+1,2])]
        count[i1, i2] <- 1 + count[i1, i2]
      }
    }
    return(count)
  }

  # Count across all edges.
  count <- emptycount
  for (edge in h$history)
    count <- count + count.edge(edge)
  # (decided against lapply; memory, and do.call frustration)

  # Cosmetics
  if (flatten) {
    count <- flatten.transition.matrix(count, sn)
  } else {
    if (is.null(sn))
        sn <- h$states
    dimnames(count) <- list(paste("fr", sn, sep="."), 
                            paste("to", sn, sep="."))
  }

  return(count)
}

## Convert a "from-to" matrix to a vector.
##   Diagonal elements are skipped.
##   Other elements are ordered nicely and named if desired..
## m = matrix; sn = state names
flatten.transition.matrix <- function(m, sn=NULL) {
  k <- nrow(m)
  stopifnot(k == ncol(m))
  kseq <- seq_len(k)
  idx.q <- cbind(rep(kseq, each=k-1), 
                 unlist(lapply(kseq, function(i) (kseq)[-i])))
  ans <- m[idx.q]
  if (!is.null(sn) & length(sn) == k)
    names(ans) <- apply(idx.q, 1, function(x) paste(sn[x], collapse=""))

  return(ans)
}

## Given a history object (h, from simulation or stochastic mapping), 
## report how much time is spent in each state.
## Not speed optimized.  But tested.
## sn = optional names for states (to over-ride the numbers in h)
summarize.durations <- function(h, phy, sn = NULL)
{
  if (class(h) != "history")
    stop("Need to provide a history object to summarize.durations()")

  edge.pieces <- list()
  for (i in seq_along(h$history))
  {
    edge <- h$history[[i]]
    # append a column with the length of each branch piece:
    edge <- cbind(edge, c(edge[-1,1], phy$edge.length[i]) - edge[,1])
    edge.pieces[[i]] <- edge
  }

  edge.pieces <- data.frame(do.call(rbind, edge.pieces)[,-1])
  names(edge.pieces) <- c("state", "length")

  ans <- sapply(h$states, function(i) 
                sum(subset(edge.pieces, state == i)$length))
  if (is.null(sn))
  {
      names(ans) <- h$states
  } else {
      names(ans) <- sn
  }

  return(ans)
}

## Given a history object (h, from simulation or stochastic mapping), 
## report the root state.
get.root.state <- function(h, sn = NULL) {
  if (class(h) != "history")
    stop("Need to provide a history object to summarize.transitions()")

  if (is.null(sn))
  {
    ans <- h$node.state[1]
  } else {
    ans <- sn[h$node.state[1]]
  }

  return(ans)
}

#--------------------------------------------------
# Plotting functions
#--------------------------------------------------

get.ci <- function(dat, varname, p=0.9, by.tree=F)
{
    if (by.tree)
    {
        cis <- by(dat, dat$tree, 
                  function(x) diversitree:::hdr(x[[varname]], p=p))
        cis <- data.frame(cbind(unlist(lapply(cis, function(x) x[1])), 
                     unlist(lapply(cis, function(x) x[2]))))
    } else {
        cis <- diversitree:::hdr(dat[[varname]], p=p)
    }
    names(cis) <- paste(p, c("Low", "Hi"), sep="")
    return(cis)
}

plot.by.tree <- function(dat, q.names, p=0.9, errcol="#555555", hline=NULL)
{
    ntrees <- length(levels(dat$tree))
    for (q.name in q.names)
    {
        # Each tree
        cis <- get.ci(dat, q.name, p, by.tree=T)
        mds <- by(dat, dat$tree, function(x) median(x[[q.name]]))
        errbar(seq(ntrees), mds, cis[,2], cis[,1], ylab="rate estimate",
               xlab="tree",  xlim=c(0, ntrees), errbar.col=errcol, col=errcol)
        title(q.name, font.main=3)
        # Overall
        cis <- get.ci(dat, q.name, p)
        errbar(0, median(dat[[q.name]]), cis[2], cis[1], add=T, lwd=4, cex=1.5)
        # Reference line
        if (!is.null(hline))
            abline(h = hline, lty=2)
    }
}

### These are from eeg_plots.R, based on diversitree ###

add.profile.outline <- function(h, col, vertical=FALSE, lwd=1, lty=1) {
  dx <- diff(h$mids[1:2])
  if ( vertical )
    lines(h, freq=FALSE, col=col, lwd=lwd, lty=lty)
  else
    with(h, lines(mids - dx/2, counts, type="s", col=col, lwd=lwd, lty=lty))
}

add.alpha <- function(col, alpha=.5) {
  tmp <- col2rgb(col)/255
  rgb(tmp[1,], tmp[2,], tmp[3,], alpha=alpha)
}

add.profile.shading <- function(h, col, ci=c(-Inf, Inf)) {
    dx <- diff(h$mids[1:2])
    i <- which(with(h, mids > ci[1] & mids < ci[2]))
    with(h, polygon(rep(c(mids[c(i, i[length(i)]+1)] - dx/2), each=2),
                    c(0, rep(counts[i], each=2), 0),
                    col=col, border=NA))
    # current diversitree has:
    # dx <- diff(h$mids[1:2])
    # xx <- c(h$mids - dx/2, h$mids[length(h$mids)] + dx/2)
    # i <- which(xx > ci[1] & xx < ci[2])
    # xs <- rep(c(ci[1], xx[i], ci[2]), each = 2)
    # j <- if (length(i) > 1) max(1, min(i) - 1) else 1
    # ys <- c(0, rep(h$density[c(j, i)], each = 2), 0)
    # polygon(xs, ys, col = col, border = NA)
}

hist.overlay <- function(hists, xlim, ylim, col.line=rep("black",999), 
                      col.fill=NULL, lty=rep(1,999), lwd=rep(3,999), cex=1.2)
{
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)

    for ( i in seq_along(hists) )
    {
        counts <- hists[[i]]$counts
        nbin <- length(counts)
        if (counts[nbin] != 0)
        {
            hists[[i]]$counts <- c(counts, 0)
            mids <- hists[[i]]$mids
            hists[[i]]$mids <- c(mids, mids[nbin]*2 - mids[nbin-1])
        }
    }


    if (!is.null(col.fill))
    {
        for ( i in seq_along(hists) )
            add.profile.shading(hists[[i]], col.fill[i])
    }
    for ( i in seq_along(hists) )
        add.profile.outline(hists[[i]], col.line[i], lwd=lwd[i], lty=lty[i])

    axis(1, cex.axis=cex)
    axis(2, cex.axis=cex)
    box()
}

hists.plot <- function(dat, col.line, legend.pos=NULL)
{
    xlim <- c(min(dat), max(dat))
    breaks <- seq(from=xlim[1], to=xlim[2]+1, by=1) - 0.5
    hists <- apply(dat, 2, hist, breaks=breaks, plot=F)
    xlim <- range(breaks)
    ylim <- c(0, max(unlist(lapply(hists, function(x) x$counts))))
    col.fill <- add.alpha(col.line)
    hist.overlay(hists, xlim, ylim, col.line=col.line, col.fill=col.fill)
    if (!is.null(legend.pos))
        legend(legend.pos, names(hists), fill=col.fill, border=col.line, bty="n")
}
