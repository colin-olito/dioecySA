#####################################################
#  Sexually antagonistic polymorphism and dioecy
#
#  R code to run deterministic simulations of the
#  genotypic frequency recursions for all four models
#  via invasion of recessive sterility alleles
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
source('R/functions-figures.R')
source('R/functions-recSim-Gyno-Rec.R')

######################################
#  Additive effects (hf = hm = 0.5)

#  Linear C-delta relation
	# dStar = 0.8

	# sm = 0.1
GynR0  <-	recursionFwdSimLoop(gen = 500000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, resolution = 0.1,
		                kMult = c(1.1, 0.975, 0.95, 0.925, 0.90), r.vals = c(0.0), threshold = 1e-7)

rm(list=ls())
#####################
##  Dependencies
source('R/functions-analyses.R')
source('R/functions-figures.R')
source('R/functions-recSim-Gyno-Rec.R')

######################################
#  Additive effects (hf = hm = 0.5)

#  Linear C-delta relation
	# dStar = 0.8

	# sm = 0.1
GynR005  <-		recursionFwdSimLoop(gen = 500000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, resolution = 0.1,
		                kMult = c(1.1, 0.975, 0.95, 0.925, 0.90), r.vals = c(0.005), threshold = 1e-7)



rm(list=ls())
#####################
##  Dependencies
source('R/functions-analyses.R')
source('R/functions-figures.R')
source('R/functions-recSim-Andro-Rec.R')

######################################
#  Additive effects (hf = hm = 0.5)

#  Linear C-delta relation
	# dStar = 0.8

	# sm = 0.1
AndR0  <-		recursionFwdSimLoop(gen = 500000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, resolution = 0.1,
		                kMult = c(1.1, 0.975, 0.95, 0.925, 0.90), r.vals = c(0.0), threshold = 1e-7)



rm(list=ls())
#####################
##  Dependencies
source('R/functions-analyses.R')
source('R/functions-figures.R')
source('R/functions-recSim-Andro-Rec.R')

######################################
#  Additive effects (hf = hm = 0.5)

#  Linear C-delta relation
	# dStar = 0.8

	# sm = 0.1
AndR005  <- 	recursionFwdSimLoop(gen = 500000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, resolution = 0.1,
		                kMult = c(1.1, 0.975, 0.95, 0.925, 0.90), r.vals = c(0.005), threshold = 1e-7)
