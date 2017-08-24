#--------------------------------------------------
# Source after common4.R, to over-write functions for bisse.
#--------------------------------------------------

rename.rates <- function(lik, sn)
{
    # Rename the transition rates, to make comparisons across genera easier
    an <- argnames(lik)
    q.old <- an[grep("q", an)]
    q.new <- c(paste("q", sn[1], sn[2], sep=""), 
               paste("q", sn[2], sn[1], sep=""))

    constr <- as.list(paste(q.old, q.new, sep=" ~ "))
    lik.ans <- constrain(lik, formulae=lapply(constr, formula),
                         extra=q.new)
    return(lik.ans)
}
