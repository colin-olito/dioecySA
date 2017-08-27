#  Functions to generate Figures for: 
#    
#  Title:  Sexually antagonistic polymorphism and the 
#          evolution of dimorphic sexual systems in 
#          hermaphrodites
#
#  Author: Colin Olito
#
#
#  NOTES: Run this file, either from terminal using Rscript,
#		  or interactively in R. This should create all the 
#		  figures needed to correctly compile the mansucript
#		  LaTeX file.  
#          

rm(list=ls())

library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)

#source('paths.R')
source('R/functions-analyses.R')
source('R/functions-figures.R')

###############
# PAPER FIGURES
###############

toPdf(Fig.1(), figPath(name='Fig1.pdf'), width=5, height=7.75)
embed_fonts(figPath(name='Fig1.pdf'))

toPdf(Fig2Alt(dfGyn = "./output/data/simResults/gyn-recess_dStar0.8_a1_sm0.1_add.csv",
              dfAnd = "./output/data/simResults/and-recess_dStar0.8_a1_sm0.1_add.csv"), 
              figPath(name='Fig2Alt.pdf'), width=7, height=4)
embed_fonts(figPath(name='Fig2Alt.pdf'))

toPdf(Fig.2(), figPath(name='Fig2.pdf'), width=7, height=7)
embed_fonts(figPath(name='Fig2.pdf'))

toPdf(Fig3Alt(), 
            figPath(name='Fig3Alt.pdf'), width=7, height=7)
embed_fonts(figPath(name='Fig3Alt.pdf'))


########################
# Supplementary Figures
########################

## Figs. S1-S6 -- Funnel plots for invasion into initially polymorphic populations
	# Gynodioecy
	  # Obligate Outcrossing
	  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-wksel-ObOut-Add-EQInv.csv", wkSel = TRUE), 
	  	          figPath(name='FigS1-Gyno-obOut-funnel.pdf'), width=9, height=10)
	  embed_fonts(figPath(name='FigS1-Gyno-obOut-funnel.pdf'))
	
	  # Low selfing, High Inbreeding Depression
	  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-wksel-partSelf-C25-delta80-Add-EQInv.csv", wkSel = TRUE), 
			          figPath(name='FigS2-Gyno-C25-d80-funnel.pdf'), width=9, height=10)
	  embed_fonts(figPath(name='FigS2-Gyno-C25-d80-funnel.pdf'))
	
	  # High selfing, Low Inbreeding Depression
	  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-wksel-partSelf-C75-delta20-Add-EQInv.csv", wkSel=TRUE), 
			          figPath(name='FigS3-Gyno-C75-d20-funnel.pdf'), width=9, height=10)
	  embed_fonts(figPath(name='FigS3-Gyno-C75-d20-funnel.pdf'))

	# Androdioecy
	  # Obligate Outcrossing
	  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/And-wksel-ObOut-Add-EQInv.csv", wkSel=TRUE), 
	  	          figPath(name='FigS4-Andro-obOut-funnel.pdf'), width=9, height=10)
	  embed_fonts(figPath(name='FigS4-Andro-obOut-funnel.pdf'))
	
	  # Low selfing, High Inbreeding Depression
	  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/And-wksel-partSelf-C25-delta80-Add-EQInv.csv", wkSel = TRUE), 
			          figPath(name='FigS5-Andro-C25-d80-funnel.pdf'), width=9, height=10)
	  embed_fonts(figPath(name='FigS5-Andro-C25-d80-funnel.pdf'))
	
	  # High selfing, Low Inbreeding Depression
	  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/And-wksel2-partSelf-C75-delta20-Add-EQInv.csv", wkSel=TRUE), 
			          figPath(name='FigS6-Andro-C75-d20-funnel.pdf'), width=9, height=10)
	  embed_fonts(figPath(name='FigS6-Andro-C75-d20-funnel.pdf'))


## Fig. S7 -- Equilibrium frequency plots showing effects of non-linear C*delta function and stronger selfing (from deterministic simulations)
toPdf(compareCDelta(), 
            figPath(name='FigS7-compareCdelta.pdf'), width=7, height=4)
embed_fonts(figPath(name='FigS7-compareCdelta.pdf'))

######################
# Exploratory Figures
######################

######################################################
##  Gynodioecy EQ invasion analyses exploratory plots 
  # Obligate Outcrossing
  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-ObOut-Add-EQInv.csv"), 
  	                    figPath(name='Gyn-ObOut-Add-EQInv.pdf'), width=8, height=10)
  embed_fonts(figPath(name='Gyn-ObOut-Add-EQInv.pdf'))


  # Partial Selfing

  # Low selfing, High Inbreeding Depression
  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-partSelf-C25-delta80-strgSel-Add-EQInv.csv", wkSel = FALSE), 
		  figPath(name='Gyn-PartSelf-Add-EQInv-C0.25-delta0.8.pdf'), width=8, height=10)
  embed_fonts(figPath(name='Gyn-PartSelf-Add-EQInv-C0.25-delta0.8.pdf'))

  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-wksel-partSelf-C25-delta80-strgSel-Add-EQInv.csv", wkSel = TRUE), 
		  figPath(name='Gyn-wkSel-PartSelf-Add-EQInv-C0.25-delta0.8.pdf'), width=9, height=10)
  embed_fonts(figPath(name='Gyn-wkSel-PartSelf-Add-EQInv-C0.25-delta0.8.pdf'))

  # High selfing, Low Inbreeding Depression
  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-partSelf-C75-delta20-strgSel-Add-EQInv.csv", wkSel=FALSE), 
		  figPath(name='Gyn-PartSelf-Add-EQInv-C0.75-delta0.2.pdf'), width=8, height=10)
  embed_fonts(figPath(name='Gyn-PartSelf-Add-EQInv-C0.75-delta0.2.pdf'))

  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-wksel-partSelf-C75-delta20-strgSel-Add-EQInv.csv", wkSel=TRUE), 
		  figPath(name='Gyn-wkSel-PartSelf-Add-EQInv-C0.75-delta0.2.pdf'), width=9, height=10)
  embed_fonts(figPath(name='Gyn-wkSel-PartSelf-Add-EQInv-C0.75-delta0.2.pdf'))


#######################################################
##  Androdioecy EQ invasion analyses exploratory plots 
  # Obligate Outcrossing
  toPdf(EQInv.ObOut.Add(df="./output/data/EQInvAnalyses/And-ObOut-Add-EQInv.csv"), 
  	                    figPath(name='And-ObOut-Add-EQInv.pdf'), width=6, height=10)
  embed_fonts(figPath(name='Andro-ObOut-Add-EQInv.pdf'))

  # Partial Selfing

  # Low selfing, High Inbreeding Depression
  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-partSelf-C25-delta80-strgSel-Add-EQInv.csv"), 
		  figPath(name='Andro-PartSelf-Add-EQInv-C0.25-delta0.8.pdf'), width=8, height=10)
  embed_fonts(figPath(name='Andro-PartSelf-Add-EQInv-C0.25-delta0.8.pdf'))

  # High selfing, Low Inbreeding Depression
  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-partSelf-C75-delta20-strgSel-Add-EQInv.csv"), 
		  figPath(name='Andro-PartSelf-Add-EQInv-C0.75-delta0.2.pdf'), width=8, height=10)
  embed_fonts(figPath(name='Andro-PartSelf-Add-EQInv-C0.75-delta0.2.pdf'))


############################################################
##  Gynodioecy dominant simulation plots
source('R/functions-recSim-Gyno-Dom.R')

# Additive fitness effects
  # sm = 0.4
  toPdf(gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.4_add.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.4_add.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.4_add.pdf'))

  # sm = 0.1
  toPdf(gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.1_add.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.1_add.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.1_add.pdf'))


# Dominance Reversal
  # sm = 0.4
  toPdf(gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.4_domRev.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.4_domRev.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.4_domRev.pdf'))

  # sm = 0.1
  toPdf(gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.1_domRev.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.1_domRev.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.1_domRev.pdf'))


############################################################
##  Androdioecy dominant simulation plots
source('R/functions-recSim-Gyno-Dom.R')

# Additive fitness effects
  # sm = 0.4
  toPdf(andDomRecPlots(df = "./output/data/simResults/and-dom_dStar0.8_a1_sm0.4_add.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.4_add.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.4_add.pdf'))

# Dominance Reversal
  # sm = 0.4
  toPdf(andDomRecPlots(df = "./output/data/simResults/and-dom_dStar0.8_a1_sm0.4_domRev.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.4_domRev.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.4_domRev.pdf'))


############################################################
##  Gynodioecy recessive simulation plots
source('R/functions-recSim-Gyno-Rec.R')

# Additive fitness effects
  # sm = 0.1
  toPdf(gynRecRecPlots(df = "./output/data/simResults/gyn-recess_dStar0.8_a1_sm0.1_add.csv"), 
  	                    figPath(name='Gyn-rec-Recursion-dStar0.8_a1_sm0.1_add.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-rec-Recursion-dStar0.8_a1_sm0.1_add.pdf'))

gynRecRecPlots(df = "./output/data/simResults/gyn-recess_dStar0.8_a1_sm0.1_add.csv")
gynRecRecPlots(df = "./output/data/simResults/gyn-recess_dStar0.8_a1_sm0.1_domRev.csv")
gynRecRecPlots(df = "./output/data/simResults/gyn-recess_dStar0.8_a0.2_sm0.1_add.csv")
gynRecRecPlots(df = "./output/data/simResults/gyn-recess_dStar0.8_a0.2_sm0.1_domRev.csv")
############################################################
##  Androdioecy recessive simulation plots
source('R/functions-recSim-Andro-Rec.R')

# Additive fitness effects
  # sm = 0.1
  toPdf(andRecRecPlots(df = "./output/data/simResults/and-recess_dStar0.8_a1_sm0.1_add.csv"), 
  	                    figPath(name='And-rec-Recursion-dStar0.8_a1_sm0.1_add.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='And-rec-Recursion-dStar0.8_a1_sm0.1_add.pdf'))


andRecRecPlots(df = "./output/data/simResults/and-recess_dStar0.8_a1_sm0.1_add.csv")
andRecRecPlots(df = "./output/data/simResults/and-recess_dStar0.8_a1_sm0.1_domRev.csv")
andRecRecPlots(df = "./output/data/simResults/and-recess_dStar0.8_a0.2_sm0.1_add.csv")
andRecRecPlots(df = "./output/data/simResults/and-recess_dStar0.8_a0.2_sm0.1_domRev.csv")
