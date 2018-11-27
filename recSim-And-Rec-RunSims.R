#####################################################
#  Sexually antagonistic polymorphism and dioecy
#
#  R code to run deterministic simulations of the
#  genotypic frequency recursions for the model of 
#  ANDRODIOECY via the invasion of a Completely 
#  RECESSIVE FEMALE sterility allele
#
#  Appendix XXX:
#  
#  Title: Sexually antagonistic polymorphism and the
#         evolution of dimorphic sexual systems in 
#         hermaphrodites
#
#  
#  Author: Colin Olito
#
#  NOTES:  
#          

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

			# all sims (TAKES A VERY LONG TIME)
			# k = kHat*1.1
#				recursionFwdSimLoop.oneLev(gen = 500000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, resolution = 0.1,
#										   kMult = c(1.1, 0.975, 0.95, 0.925, 0.90), r.vals = c(0.0, 0.005, 0.01, 0.05), threshold = 1e-7)

			# Single values for kMult and r
			# k = kHat*1.1
				recursionFwdSimLoop.oneLev(gen = 500000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, resolution = 0.1,
										   kMult = c(1.1), r.vals = c(0.0, 0.005, 0.01, 0.05), threshold = 1e-7)
			# k = kHat*1.0
				recursionFwdSimLoop.oneLev(gen = 500000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, resolution = 0.1,
										   kMult = c(1.001), r.vals = c(0.0, 0.005, 0.01, 0.05), threshold = 1e-7)
			# r = 0
				recursionFwdSimLoop.oneLev(gen = 500000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, resolution = 0.1,
										   kMult = c(1.1, 0.975, 0.95, 0.925, 0.90), r.vals = c(0.0), threshold = 1e-7)
			# r = 0.005
				recursionFwdSimLoop.oneLev(gen = 500000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, resolution = 0.1,
										   kMult = c(1.1, 0.975, 0.95, 0.925, 0.90), r.vals = c(0.005), threshold = 1e-7)


#  Concave C-delta relation
	# dStar = 0.8

	# sm = 0.1
	recursionFwdSimLoop(gen = 25000, dStar = 0.8, a = 0.2, b = 0.5, sm = 0.1, hf = 0.5, hm = 0.5, resolution = 0.1,
		                kMult = c(1.1, 0.975, 0.95, 0.925, 0.90), r.vals = c(0.0, 0.005, 0.01, 0.05), threshold = 1e-7)



######################################
#  Dominance Reversal (hf = hm = 0.25)

#  Linear C-delta relation
	# dStar = 0.8

	# sm = 0.1
	recursionFwdSimLoop(gen = 25000, dStar = 0.8, a = 1, b = 0.5, sm = 0.1, hf = 0.25, hm = 0.25, resolution = 0.1, 
		                kMult = c(1.1, 0.975, 0.95, 0.925, 0.90), r.vals = c(0.0, 0.005, 0.01, 0.05), threshold = 1e-7)


#  Concave C-delta relation
	# dStar = 0.8

	# sm = 0.1
	recursionFwdSimLoop(gen = 25000, dStar = 0.8, a = 0.2, b = 0.5, sm = 0.1, hf = 0.25, hm = 0.25, resolution = 0.1, 
		                kMult = c(1.1, 0.975, 0.95, 0.925, 0.90), r.vals = c(0.0, 0.005, 0.01, 0.05), threshold = 1e-7)
