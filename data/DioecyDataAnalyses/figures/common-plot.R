# Threshold for statistical significance
p <- 0.95

# Colors for forward-backward

myblue <- rgb(0, 0.445, 0.695)      # 0071B1
myorange <- rgb(0.832, 0.367, 0)    # D45E00

# take colors above, set value=100 and saturation=42
# gimp or colorizer.org
mylightblue <- rgb(0.576, 0.847, 1)
mylightorange <- rgb(1, 0.765, 0.576)

col.dark <- c(myblue, myorange)
col.light <- c(mylightblue, mylightorange)
rm(myblue, myorange, mylightblue, mylightorange)

mygray <- rgb(0.7, 0.7, 0.7)

# Colors for four states

col.state <- c(H = "gold1", M = "green4", G = "blue2", D = "red4", O = "gray")

# Colors for root state

rootH <- c("Lycium", "Thalictrum", "Cordia")
rootD <- c("Bursera", "Dodonaea")
rootM <- c("Croton", "Begonia")
# no rootG

col.root <- c("H" = rgb(1, 0.933, 0.580),
              "M" = rgb(0.486, 0.839, 0.486),
              "G" = rgb(0.541, 0.541, 0.933),
              "D" = rgb(0.8, 0.463, 0.463)
              )

# Large genera
genera.large <- read.table("genera-large.dat", as.is=T)[,1]

# Interleave elements of two vectors, without repeating.
# https://stat.ethz.ch/pipermail/r-help/2010-March/233149.html
riffle <- function(a, b)
{
    n <- min(length(a), length(b))
    p1 <- as.vector(rbind(a[1:n], b[1:n]))
    p2 <- c(a[-(1:n)], b[-(1:n)])
    c(p1, p2)
}
