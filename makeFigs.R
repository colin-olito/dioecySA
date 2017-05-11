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

toPdf(Fig.2(), figPath(name='Fig2.pdf'), width=7, height=7)
embed_fonts(figPath(name='Fig2.pdf'))

#toPdf(Fig.2wk(), figPath(name='Fig2wk.pdf'), width=7, height=7)
#embed_fonts(figPath(name='Fig2wk.pdf'))

#toPdf(recSimFig_add(), figPath(name='recSimFig_add.pdf'), width=7, height=7)
#embed_fonts(figPath(name='recSimFig_add.pdf'))

#toPdf(recSimFig_domRev(), figPath(name='recSimFig_domRev.pdf'), width=7, height=7)
#embed_fonts(figPath(name='recSimFig_domRev.pdf'))


########################
# Supplementary Figures
########################

# Figs. S1-S6 -- Funnel plots for invasion into initially polymorphic populations
	# Gynodioecy
	  # Obligate Outcrossing
	  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-wksel-ObOut-Add-EQInv.csv", wkSel = TRUE), 
	  	          figPath(name='FigS1-Gyno-obOut-funnel.pdf'), width=8, height=10)
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
	  	          figPath(name='FigS4-Andro-obOut-funnel.pdf'), width=8, height=10)
	  embed_fonts(figPath(name='FigS4-Andro-obOut-funnel.pdf'))
	
	  # Low selfing, High Inbreeding Depression
	  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/And-wksel-partSelf-C25-delta80-Add-EQInv.csv", wkSel = TRUE), 
			      figPath(name='FigS5-Andro-C25-d80-funnel.pdf'), width=9, height=10)
	  embed_fonts(figPath(name='FigS5-Andro-C25-d80-funnel.pdf'))
	
	  # High selfing, Low Inbreeding Depression
	  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-wksel-partSelf-C75-delta20-Add-EQInv.csv", wkSel=TRUE), 
			      figPath(name='FigS6-Andro-C75-d20-funnel.pdf'), width=9, height=10)
	  embed_fonts(figPath(name='FigS6-Andro-C75-d20-funnel.pdf'))



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
#  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-partSelf-C25-delta20-strgSel-Add-EQInv.csv"), 
#     	  figPath(name='Gyn-PartSelf-Add-EQInv-C0.25-delta0.2.pdf'), width=8, height=10)
#  embed_fonts(figPath(name='Gyn-PartSelf-Add-EQInv-C0.25-delta0.2.pdf'))

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

#  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-partSelf-C75-delta80-strgSel-Add-EQInv.csv"), 
#		  figPath(name='Gyn-PartSelf-Add-EQInv-C0.75-delta0.8.pdf'), width=8, height=10)
#  embed_fonts(figPath(name='Gyn-PartSelf-Add-EQInv-C0.75-delta0.8.pdf'))

#######################################################
##  Androdioecy EQ invasion analyses exploratory plots 
  # Obligate Outcrossing
  toPdf(EQInv.ObOut.Add(df="./output/data/EQInvAnalyses/And-ObOut-Add-EQInv.csv"), 
  	                    figPath(name='And-ObOut-Add-EQInv.pdf'), width=6, height=10)
  embed_fonts(figPath(name='Andro-ObOut-Add-EQInv.pdf'))

  # Partial Selfing
#  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-partSelf-C25-delta20-strgSel-Add-EQInv.csv"), 
#     	  figPath(name='Andro-PartSelf-Add-EQInv-C0.25-delta0.2.pdf'), width=8, height=10)
#  embed_fonts(figPath(name='Andro-PartSelf-Add-EQInv-C0.25-delta0.2.pdf'))

  # Low selfing, High Inbreeding Depression
  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-partSelf-C25-delta80-strgSel-Add-EQInv.csv"), 
		  figPath(name='Andro-PartSelf-Add-EQInv-C0.25-delta0.8.pdf'), width=8, height=10)
  embed_fonts(figPath(name='Andro-PartSelf-Add-EQInv-C0.25-delta0.8.pdf'))

  # High selfing, Low Inbreeding Depression
  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-partSelf-C75-delta20-strgSel-Add-EQInv.csv"), 
		  figPath(name='Andro-PartSelf-Add-EQInv-C0.75-delta0.2.pdf'), width=8, height=10)
  embed_fonts(figPath(name='Andro-PartSelf-Add-EQInv-C0.75-delta0.2.pdf'))

#  toPdf(EQInv.Add(df="./output/data/EQInvAnalyses/Gyn-partSelf-C75-delta80-strgSel-Add-EQInv.csv"), 
#		  figPath(name='Andro-PartSelf-Add-EQInv-C0.75-delta0.8.pdf'), width=8, height=10)
#  embed_fonts(figPath(name='Andro-PartSelf-Add-EQInv-C0.75-delta0.8.pdf'))


############################################################
##  Gynodioecy dominant simulation plots
source('R/functions-recSim-Gyno-Dom.R')

# Additive fitness effects
  # sm = 0.4
  toPdf(gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.4_add.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_add.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_add.pdf'))

  # sm = 0.3
  toPdf(gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.3_add.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_add.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_add.pdf'))

  # sm = 0.1
  toPdf(gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.1_add.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.1_add.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.1_add.pdf'))


# Dominance Reversal
  # sm = 0.4
  toPdf(gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.5_domRev.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_domRev.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_domRev.pdf'))



gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.4_add.csv")
gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a0.2_sm0.4_add.csv")
gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.4_domRev.csv")



############################################################
##  Androdioecy dominant simulation plots
source('R/functions-recSim-Gyno-Dom.R')

# Additive fitness effects
  # sm = 0.4
  toPdf(andDomRecPlots(df = "./output/data/simResults/and-dom_dStar0.8_a1_sm0.4_add.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_add.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_add.pdf'))

# Dominance Reversal
  # sm = 0.4
  toPdf(andDomRecPlots(df = "./output/data/simResults/and-dom_dStar0.8_a1_sm0.5_domRev.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_domRev.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_domRev.pdf'))



andDomRecPlots(df = "./output/data/simResults/and-dom_dStar0.8_a1_sm0.4_add.csv")


andDomRecPlots(df = "./output/data/simResults/and-dom_dStar0.8_a1_sm0.4_domRev.csv")




############################################################
##  Gynodioecy dominant simulation plots
source('R/functions-recSim-Gyno-Rec.R')

# Additive fitness effects
  # sm = 0.4
  toPdf(gynDomRecPlots(df = "./output/data/simResults/gyn-dom_dStar0.8_a1_sm0.4_add.csv"), 
  	                    figPath(name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_add.pdf'), width=6, height=6)
  embed_fonts(figPath(          name='Gyn-Dom-Recursion-dStar0.8_a1_sm0.5_add.pdf'))

gynRecRecPlots(df = "./output/data/simResults/gyn-recess_dStar0.8_a1_sm0.4_add.csv")
