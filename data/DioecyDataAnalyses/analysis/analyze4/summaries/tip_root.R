### Get the available genera

tmp <- system("ls -d ../[A-Z]*", intern=T)
genera <- sapply(tmp, function(x) rev(strsplit(x, "/")[[1]])[1])
names(genera) <- NULL
# Should be able to do this with dir(...pattern), but oh well.

#--------------------------------------------------
# Tip state frequencies
#--------------------------------------------------

### Construct the table framework

state.names <- c("Herm", "xMono", "xDio", "Dio")
freq.table <- matrix(0, nrow=length(genera), ncol=4,
                     dimnames=list(genera, state.names))

for (genus in genera)
{
    dat <- read.csv(paste("../../../data/", genus, ".csv", sep=""), 
                    as.is=T)$SexSyst4
    states <- unlist(strsplit(dat, split="|", fixed=T))
    ans <- table(states)
    freq.table[genus, names(ans)] <- ans
}

colnames(freq.table) <- c("H", "M", "G", "D")
write.csv(freq.table, file="tip_states.csv")

#--------------------------------------------------
# Construct the table of root state probabilities
# Construct the table of time spent in each state
#--------------------------------------------------

X.all <- c("D", "G", "M", "H")

root.table <- matrix(NA, nrow=length(genera), ncol=length(X.all),
                       dimnames=list(genera, X.all))
time.table <- matrix(NA, nrow=length(genera), ncol=length(X.all),
                       dimnames=list(genera, X.all))

for (genus in genera)
{
    here.dir <- setwd(paste("../", genus, sep=""))

    infiles <- Sys.glob("maps_tree*.csv")
    stopifnot(length(infiles) > 0)

    for (i in seq_along(infiles))
    {
        dat1 <- read.csv(infiles[i], as.is=T)
        junk <- ifelse(i == 1, dat <- dat1, dat <- rbind(dat, dat1))
    }

    ans <- table(dat$root) / nrow(dat)
    root.table[genus,] <- ans[colnames(root.table)]

    X.have <- X.all[X.all %in% names(dat)]
    ans <- colSums(dat[,X.have]) / sum(dat[,X.have])
    time.table[genus,] <- ans[colnames(time.table)]

    rm(dat1, junk, ans)
    setwd(here.dir)
    message(paste("done with", genus))
}

write.csv(root.table, file="root_states.csv")
write.csv(time.table, file="state_durations.csv")
