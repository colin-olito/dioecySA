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
source('R/functions-figures.R')
source('R/functions-recSim-Gyno-Dom.R')

######################################
#  Additive effects (hf = hm = 0.5)

#  Linear C-delta relation
	# dStar = 0.8

	# sm = 0.4
	recursionFwdSimLoop(gen = 10000, dStar = 0.8, a = 1, b = 0.5, sm = 0.4, hf = 0.5, hm = 0.5, 
		                kMult = c(1.1, 0.95, 0.90), r.vals = c(0.0, 0.05, 0.1), threshold = 1e-7)

	# sm = 0.3
	recursionFwdSimLoop(gen = 10000, dStar = 0.8, a = 1, b = 0.5, sm = 0.4, hf = 0.5, hm = 0.5, 
		                kMult = c(1.1, 0.95, 0.90), r.vals = c(0.0, 0.05, 0.1), threshold = 1e-7)

	# sm = 0.1
	recursionFwdSimLoop(gen = 15000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, 
		                kMult = c(1.1, 0.99, 0.98), r.vals = c(0.0, 0.005), threshold = 1e-6)


#  Concave C-delta relation
	# dStar = 0.8

	# sm = 0.4
	recursionFwdSimLoop(gen = 10000, dStar = 0.8, a = 0.2, b = 0.5, sm = 0.4, hf = 0.5, hm = 0.5, 
		                kMult = c(1.1, 0.95, 0.90), r.vals = c(0.0, 0.05, 0.1), threshold = 1e-7)

	# sm = 0.3
	recursionFwdSimLoop(gen = 10000, dStar = 0.8, a = 0.2, b = 0.5, sm = 0.4, hf = 0.5, hm = 0.5, 
		                kMult = c(1.1, 0.95, 0.90), r.vals = c(0.0, 0.05, 0.1), threshold = 1e-7)

	# sm = 0.1
	recursionFwdSimLoop(gen = 15000, dStar = 0.8, a = 0.2, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, 
		                kMult = c(1.1, 0.99, 0.98), r.vals = c(0.0, 0.005), threshold = 1e-6)



######################################
#  Dominance Reversal (hf = hm = 0.25)


#  Linear C-delta relation
	# dStar = 0.8

	# sm = 0.4
	recursionFwdSimLoop(gen = 10000, dStar = 0.8, a = 1, b = 0.5, sm = 0.4, hf = 0.25, hm = 0.25, 
		                kMult = c(1.1, 0.95, 0.90), r.vals = c(0.0, 0.05, 0.1), threshold = 1e-7)

	# sm = 0.3
	recursionFwdSimLoop(gen = 10000, dStar = 0.8, a = 1, b = 0.5, sm = 0.4, hf = 0.25, hm = 0.25, 
		                kMult = c(1.1, 0.95, 0.90), r.vals = c(0.0, 0.05, 0.1), threshold = 1e-7)

	# sm = 0.1
	recursionFwdSimLoop(gen = 15000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.25, hm = 0.25, 
		                kMult = c(1.1, 0.99, 0.98), r.vals = c(0.0, 0.005), threshold = 1e-6)


#  Concave C-delta relation
	# dStar = 0.8

	# sm = 0.4
	recursionFwdSimLoop(gen = 10000, dStar = 0.8, a = 0.2, b = 0.5, sm = 0.4, hf = 0.25, hm = 0.25, 
		                kMult = c(1.1, 0.95, 0.90), r.vals = c(0.0, 0.05, 0.1), threshold = 1e-7)

	# sm = 0.3
	recursionFwdSimLoop(gen = 10000, dStar = 0.8, a = 0.2, b = 0.5, sm = 0.4, hf = 0.25, hm = 0.25, 
		                kMult = c(1.1, 0.95, 0.90), r.vals = c(0.0, 0.05, 0.1), threshold = 1e-7)

	# sm = 0.1
	recursionFwdSimLoop(gen = 15000, dStar = 0.8, a = 0.2, b = 0.5, sm = 0.1, hf = 0.25, hm = 0.25, 
		                kMult = c(1.1, 0.99, 0.98), r.vals = c(0.0, 0.005), threshold = 1e-6)
