#####################################################
#  Reconsidering the evolution of dioecy from the
#  perspective of sexually antagonistic selection
#
#  Necessary functions for analysis of a model of
#  the evolution of gynodioecy via invasion of a
#  completely dominant male sterility allele 
#
#  Author: Colin Olito
#
#  NOTES:  
#          



###################################
##  A few basic invasion criteria


#' Invasion criteria for SA alleles in the case of obligate outcrossing 
#' 
#' @title Invasion criteria for SA alleles in the case of obligate outcrossing.
#' @param sm values for the strength of selection through male sex-function
#' where 0 <= sm <= 1.
#' @return Returns a vector of predicted sf values at which the opposite SA
#' allele is able to invade an obligately outcrossing population.
#' @seealso `recursionFwdSim`
#' @export
#' @author Colin Olito.
inv.AKid  <-  function(sm) {
	sm / (1 - sm)
}
inv.aKid  <-  function(sm) {
	sm / (1 + sm)
}