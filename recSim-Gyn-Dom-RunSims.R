#####################################################
#  Sexually antagonistic polymorphism and dioecy
#
#  R code to run deterministic simulations of the
#  genotypic frequency recursions for the model of 
#  Gynodioecy via the invasion of a Completely 
#  Dominant male sterility allele
#
#  Appendix XXX:
#  
#  Title: Sexually antagonistic polymorphism and the
#         evolution of dimorphic sexual systems in 
#         hermaphrodites
#
#  
#  Authors: Colin Olito, Tim Connallon
#
#  NOTES:  
#          



rm(list=ls())
#####################
##  Dependencies
source('R/functions-analyses.R')
source('R/functions-recSim-Gyno-Dom.R')

# set seed for random number generator
randSeed  <-  3497016


######################################
#  Strong selection (sm,sf ~ [0,1])
######################################

######################################
#  Additive effects (hf = hm = 0.5)

#  Obligate Outcrossing (C = 0)
recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0, delta = 0, 
	                       hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

#  Partial Selfing
recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0.25, delta = 0.2, 
	                       hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0.25, delta = 0.8, 
	                       hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0.75, delta = 0.2, 
	                       hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0.75, delta = 0.8, 
	                       hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

######################################
#  Dominance Reversal (hf = hm = 0.25)

#  Obligate Outcrossing (C = 0)
recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0, delta = 0, 
	                       hf = 0.25, hm = 0.25, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

#  Partial Selfing
recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0.25, delta = 0.2, 
	                       hf = 0.25, hm = 0.25, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0.25, delta = 0.8, 
	                       hf = 0.25, hm = 0.25, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0.75, delta = 0.2, 
	                       hf = 0.25, hm = 0.25, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0.75, delta = 0.8, 
	                       hf = 0.25, hm = 0.25, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)


######################################
#  Weak Selection (sm,sf ~ [0,0.1])
######################################

######################################
#  Additive effects (hf = hm = 0.5)

#  Obligate Outcrossing (C = 0)
recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,0.1), C = 0, delta = 0, 
	                       hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

#  Partial Selfing
recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,0.1), C = 0.25, delta = 0.2, 
	                       hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,0.1), C = 0.25, delta = 0.8, 
	                       hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,0.1), C = 0.75, delta = 0.2, 
	                       hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,0.1), C = 0.75, delta = 0.8, 
	                       hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

######################################
#  Dominance Reversal (hf = hm = 0.25)

#  Obligate Outcrossing (C = 0)
recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,0.1), C = 0, delta = 0, 
	                       hf = 0.25, hm = 0.25, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

#  Partial Selfing
recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,0.1), C = 0.25, delta = 0.2, 
	                       hf = 0.25, hm = 0.25, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,0.1), C = 0.25, delta = 0.8, 
	                       hf = 0.25, hm = 0.25, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,0.1), C = 0.75, delta = 0.2, 
	                       hf = 0.25, hm = 0.25, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)

recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,0.1), C = 0.75, delta = 0.8, 
	                       hf = 0.25, hm = 0.25, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                       seed = randSeed, threshold = 1e-7)
