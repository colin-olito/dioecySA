options(stringsAsFactors = FALSE)
require(diversitree, quietly=TRUE)
require(parallel, quietly=TRUE)
require(stringr, quietly=TRUE)
source("../../scripts/check_states.R") # patch for diversitree

## Collapse to four standard categories
#       Herm:  hermaphroditism  (all cosexual/perfect flowers)
#       Dio:   dioecy           (each individual is either male or female)
#       xMono: mixy-moo monoecy (all individuals have both sexes)
#       xDio:  mixy-moo dioecy  (individuals differ in their sexes)
traits4 <- c("Herm", "xMono", "xDio", "Dio")

#--------------------------------------------------
# Delete species
#--------------------------------------------------
## Provide a vector or list of species names.

drop.species <- function(x, ...) UseMethod("drop.species", x)

drop.species.default <- function(dat, dropme, quiet = FALSE)
{
    dropme <- clean.dropme(dropme, quiet=quiet)

    # If var/ssp is provided, drop only that and not the whole species.
    i <- which(dat$NameInfr %in% dropme)
    if (length(i) == 0)
        i <- which(dat$Name %in% dropme)
    if (length(i) > 0)
    {
        message(paste("Dropping", length(i), "species from trait matrix."))
        dat <- dat[-i,]
    } else {
        message("No species found to drop from trait matrix.")
    }

    return(dat)
}

drop.species.phylo <- function(phy, dropme, quiet=FALSE)
{
    dropme <- clean.dropme(dropme, quiet=quiet)

    # If the tips to drop aren't present, don't return a useless empty tree.
    if (any(dropme %in% phy$tip.label))
    {
        if (!quiet)
            message(paste("Dropping", sum(dropme %in% phy$tip.label), 
                          "species from phylogeny."))
        phy <- diversitree:::drop.tip.fixed(phy, dropme)
    } else {
        if (!quiet)
            message("No species found to drop from phylogeny.")
    }

    return(phy)
}

clean.dropme <- function(dropme, mark.infr = FALSE, quiet = FALSE)
{
    dropme <- unique(unlist(dropme))
    dropme <- dropme[dropme != ""]
    if (length(dropme) == 0 & !quiet)
        message("No species provided to drop.")

    # Note which dropme names have infraspecific modifiers.
    if (mark.infr)
    {
        if (length(dropme) == 0)
        {
            dropme <- list(infr = c(), sp = c())
        } else {
            # These are the dropme's that require infraspecific matches.
            i <- grep.full(infraspecific.modifiers, dropme)

            if (length(i) > 0)
            {
                dropme <- list(infr = dropme[i], sp = dropme[-i])
            } else {
                dropme <- list(infr = c(), sp = dropme)
            }
        }
    }

    return(dropme)
}

#--------------------------------------------------
# General utilities
#--------------------------------------------------

# Reverse the order of letters in a string
rev.str <- function(a)
{
    paste(rev(substring(a, 1:nchar(a),1:nchar(a))), collapse="")
}

# Fill with leading zeros to get an n-digit character string
zfill <- function(x, n)
{
    nc <- nchar(x)
    zeros <- paste(rep(0, n), collapse = "")
    paste(substring(zeros, nchar(x) + 1, n), substring(x, 1, nchar(x)), 
          sep = "")
}

# Search for *any* of the elements of patterns.
grep.full <- function(patterns, x, ...)
{
    patterns <- unique(patterns)
    i <- c()
    if (length(patterns) > 0) 
    {
        # (faster than lapply/unlist)
        for (p in patterns)
            i <- c(i, grep(p, x, ...))
        i <- sort(unique(i))
    }
    return(i)
}

## Preferred defaults for writing CSV files.
write.csv0 <- function(x, filename)
{
    write.csv(x, file = filename, row.names = FALSE, na = "")
}
