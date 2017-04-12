#  Functions to generate Figures for: 
#    
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

#toPdf(Fig.1(), figPath(name='Fig1.pdf'), width=7, height=7)
#embed_fonts(figPath(name='Fig1.pdf'))

#toPdf(Fig.1wk(), figPath(name='Fig1wk.pdf'), width=7, height=7)
#embed_fonts(figPath(name='Fig1wk.pdf'))

#toPdf(Fig.2(), figPath(name='Fig2.pdf'), width=7, height=7)
#embed_fonts(figPath(name='Fig2.pdf'))

#toPdf(Fig.2wk(), figPath(name='Fig2wk.pdf'), width=7, height=7)
#embed_fonts(figPath(name='Fig2wk.pdf'))

#toPdf(recSimFig_add(), figPath(name='recSimFig_add.pdf'), width=7, height=7)
#embed_fonts(figPath(name='recSimFig_add.pdf'))

#toPdf(recSimFig_domRev(), figPath(name='recSimFig_domRev.pdf'), width=7, height=7)
#embed_fonts(figPath(name='recSimFig_domRev.pdf'))



######################
# Exploratory Figures
######################

toPdf(EQInv.ObOut.Add.Gyno.Dom(), figPath(name='EQInv-ObOut-Add-Gyno-Dom.pdf'), width=6, height=10)
embed_fonts(figPath(name='EQInv-ObOut-Add-Gyno-Dom.pdf'))