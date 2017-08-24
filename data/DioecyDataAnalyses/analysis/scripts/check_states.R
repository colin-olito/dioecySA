### This modified version of "check.states" is needed to work with diversitree
### 0.9-7. It allows Mkn models to run when not all states are present.

patched.check.states <- function(tree, states, allow.unnamed=FALSE,
                         strict=FALSE, strict.vals=NULL,
                         as.integer=TRUE) {
  multicheck <- TRUE # for multistate strict checking
  if ( is.matrix(states) ) {
    ## Multistate characters (experimental).  This will not work with
    ## clade trees, but they are only interesting for BiSSE, which has
    ## NA values for multistate (even weight).
    if ( inherits(tree, "clade.tree") )
      stop("Clade trees won't work with multistate tips yet")
    n <- rowSums(states > 0)
    if ( any(n == 0) )
      stop(sprintf("No state found for taxa: %s",
                   paste(names(n)[n == 0], collapse=", ")))
    if (any(rowSums(states) == 0))
        multicheck <- FALSE

    i.mono <- which(n == 1)
    i.mult <- which(n >  1)

    tmp <- diversitree:::matrix.to.list(states)
    names(tmp) <- rownames(states)

    states.mult <- lapply(tmp[i.mult], as.numeric)

    states <- rep(NA, length(tmp))
    names(states) <- names(tmp)
    states[i.mono] <- sapply(tmp[i.mono], function(x)
                             which(x != 0))

    attr(states, "multistate") <- list(i=i.mult, states=states.mult)
  }
  
  if ( is.null(names(states)) ) {
    if ( allow.unnamed ) {
      if ( length(states) == length(tree$tip.label) ) {
        names(states) <- tree$tip.label
        warning("Assuming states are in tree$tip.label order")
      } else {
        stop(sprintf("Invalid states length (expected %d)",
                     length(tree$tip.label)))
      }
    } else {
      stop("The states vector must contain names")
    }
  }
  
  if ( !all(tree$tip.label %in% names(states)) )
    stop("Not all species have state information")

  ## When multistate characters are present, this may fail even
  ## for cases where it should not.
  ## now, multicheck helps this
  if ( !is.null(strict.vals) ) {
    if ( isTRUE(all.equal(strict.vals, 0:1)) )
      if ( is.logical(states) )
        states[] <- as.integer(states)
    
    if ( strict ) {
      if ( !isTRUE(all.equal(sort(strict.vals),
                             sort(unique(na.omit(states))))) & !multicheck)
        stop("Because strict state checking requested, all (and only) ",
             sprintf("states in %s are allowed",
                     paste(strict.vals, collapse=", ")))
    } else {
      extra <- setdiff(sort(unique(na.omit(states))), strict.vals)
      if ( length(extra) > 0 )
        stop(sprintf("Unknown states %s not allowed in states vector",
                     paste(extra, collapse=", ")))
    }
    if ( as.integer && any(!is.na(states)) )
      states <- diversitree:::check.integer(states)
  }

  if ( inherits(tree, "clade.tree") ) {
    spp.clades <- unlist(tree$clades)
    if ( !all(spp.clades %in% names(states)) )
      stop("Species in 'clades' do not have states information")
    states[union(tree$tip.label, spp.clades)]
  } else {
    ret <- states[tree$tip.label]
    ## Ugly hack...
    attr(ret, "multistate") <- attr(states, "multistate")
    ret
  }
}

assignInNamespace("check.states", patched.check.states, "diversitree")
rm(patched.check.states)
