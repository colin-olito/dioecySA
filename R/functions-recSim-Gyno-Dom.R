#####################################################
#  Reconsidering the evolution of dioecy from the
#  perspective of sexually antagonistic selection
#
#  Necessary functions for deterministic simulation
#  of genotypic frequency recursions for a model of
#  the evolution of gynodioecy via invasion of a
#  completely dominant male sterility allele 
#
#  Author: Colin Olito
#
#  NOTES:  
#          

##########################
##  Key to Fii subscripts
#  F11  =  Fii[1]
#  F12  =  Fii[2]
#  F13  =  Fii[3]
#  F14  =  Fii[4]
#  F22  =  Fii[5]
#  F23  =  Fii[6]
#  F24  =  Fii[7]
#  F33  =  Fii[8]
#  F34  =  Fii[9]
#  F44  =  Fii[10]

#######################################################
## Necessary Functions
#######################################################

# Average fitness of adults after inbreeding depression, but prior to mating
Dbar  <-  function(Fii, Gii, par.list, ...) {
	(Fii[1] + Fii[2] + Fii[3] + Fii[4] + Fii[5] + Fii[6] + Fii[7] + Fii[8] + Fii[9] + Fii[10]) + (1 - par.list$delta)*(Gii[1] + Gii[2] + Gii[3] + Gii[4] + Gii[5] + Gii[6] + Gii[7] + Gii[8] + Gii[9] + Gii[10])
}

# Frequencies of adult genotypes after inbreeding depression 

FA.11  <-  function(Fii, Gii, par.list, ...) {
	(Fii[1] + Gii[1] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.12  <-  function(Fii, Gii, par.list, ...) {
	(Fii[2] + Gii[2] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.13  <-  function(Fii, Gii, par.list, ...) {
	(Fii[3] + Gii[3] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.14  <-  function(Fii, Gii, par.list, ...) {
	(Fii[4] + Gii[4] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.22  <-  function(Fii, Gii, par.list, ...) {
	(Fii[5] + Gii[5] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.23  <-  function(Fii, Gii, par.list, ...) {
	(Fii[6] + Gii[6] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.24  <-  function(Fii, Gii, par.list, ...) {
	(Fii[7] + Gii[7] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.33  <-  function(Fii, Gii, par.list, ...) {
	(Fii[8] + Gii[8] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.34  <-  function(Fii, Gii, par.list, ...) {
	(Fii[9] + Gii[9] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.44  <-  function(Fii, Gii, par.list, ...) {
	(Fii[10] + Gii[10] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}

# Total ovules produced 
OTot  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	FA.11(Fii, Gii, par.list)*Wf.mat[1,1] + FA.12(Fii, Gii, par.list)*Wf.mat[1,2] + FA.13(Fii, Gii, par.list)*Wf.mat[1,3] + 
	FA.14(Fii, Gii, par.list)*Wf.mat[1,4] + FA.22(Fii, Gii, par.list)*Wf.mat[2,2] + FA.23(Fii, Gii, par.list)*Wf.mat[2,3] +
	FA.24(Fii, Gii, par.list)*Wf.mat[2,4] + FA.33(Fii, Gii, par.list)*Wf.mat[3,3] + FA.34(Fii, Gii, par.list)*Wf.mat[3,4] + 
	FA.44(Fii, Gii, par.list)*Wf.mat[4,4]
}

# Total ovule used for self-fertilization
OsTot  <-  function(Fii, Gii, par.list, Wf.mat) {
	if(par.list$C == 0) 0 
		else(par.list$C * (FA.11(Fii, Gii, par.list)*Wf.mat[1,1] + FA.13(Fii, Gii, par.list)*Wf.mat[1,3] + FA.33(Fii, Gii, par.list)*Wf.mat[3,3]))
}

# Total ovules used for outcrossing
OxTot  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(1 - par.list$C)*(FA.11(Fii, Gii, par.list)*Wf.mat[1,1] + FA.13(Fii, Gii, par.list)*Wf.mat[1,3] + FA.33(Fii, Gii, par.list)*Wf.mat[3,3]) + 
	(FA.12(Fii, Gii, par.list)*Wf.mat[1,2] + FA.14(Fii, Gii, par.list)*Wf.mat[1,4] + FA.22(Fii, Gii, par.list)*Wf.mat[2,2] + 
	 FA.23(Fii, Gii, par.list)*Wf.mat[2,3] + FA.24(Fii, Gii, par.list)*Wf.mat[2,4] + FA.34(Fii, Gii, par.list)*Wf.mat[3,4] + 
	 FA.44(Fii, Gii, par.list)*Wf.mat[4,4])
}

# Proportion of all ovules that are self-fertilized
S  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	OsTot(Fii, Gii, par.list, Wf.mat)/OTot(Fii, Gii, par.list, Wf.mat)
}

# Proportional contribution to selfed offspring for each genotype capable of selfing
os.11  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	if(par.list$C == 0) 0
		else(par.list$C * ((FA.11(Fii, Gii, par.list)*Wf.mat[1,1]) / OsTot(Fii, Gii, par.list, Wf.mat)))
}
os.13  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	if(par.list$C == 0) 0
		else(par.list$C * ((FA.11(Fii, Gii, par.list)*Wf.mat[1,3]) / OsTot(Fii, Gii, par.list, Wf.mat)))
}
os.33  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	if(par.list$C == 0) 0
	else(par.list$C * ((FA.11(Fii, Gii, par.list)*Wf.mat[3,3]) / OsTot(Fii, Gii, par.list, Wf.mat)))
}

# Proportional contribution to outcrossed offspring for each genotype capable of outcrossing
ox.11  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(1 - par.list$C) * ((FA.11(Fii, Gii, par.list)*Wf.mat[1,1]) / OxTot(Fii, Gii, par.list, Wf.mat))
}
ox.12  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(FA.12(Fii, Gii, par.list)*Wf.mat[1,2]) / OxTot(Fii, Gii, par.list, Wf.mat)
}
ox.13  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(1 - par.list$C) * ((FA.13(Fii, Gii, par.list)*Wf.mat[1,3]) / OxTot(Fii, Gii, par.list, Wf.mat))
}
ox.14  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(FA.14(Fii, Gii, par.list)*Wf.mat[1,4]) / OxTot(Fii, Gii, par.list, Wf.mat)
}
ox.22  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(FA.22(Fii, Gii, par.list)*Wf.mat[2,2]) / OxTot(Fii, Gii, par.list, Wf.mat)
}
ox.23  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(FA.23(Fii, Gii, par.list)*Wf.mat[2,3]) / OxTot(Fii, Gii, par.list, Wf.mat)
}
ox.24  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(FA.24(Fii, Gii, par.list)*Wf.mat[2,4]) / OxTot(Fii, Gii, par.list, Wf.mat)
}
ox.33  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(1 - par.list$C) * ((FA.33(Fii, Gii, par.list)*Wf.mat[3,3]) / OxTot(Fii, Gii, par.list, Wf.mat))
}
ox.34  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(FA.34(Fii, Gii, par.list)*Wf.mat[3,4]) / OxTot(Fii, Gii, par.list, Wf.mat)
}
ox.44  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(FA.44(Fii, Gii, par.list)*Wf.mat[4,4]) / OxTot(Fii, Gii, par.list, Wf.mat)
}

# Total pollen used for outcrossing (proportional to) the total amont of pollen 
# in the population (the amount of pollen used for selfing has negligable effect 
# on the pool of outcrossing pollen)
PxTot  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	FA.11(Fii, Gii, par.list)* Wm.mat[1,1] + FA.13(Fii, Gii, par.list)* Wm.mat[1,3] + FA.33(Fii, Gii, par.list)* Wm.mat[3,3]
}

# Frequency of pollen/sperm | Outcrossing genotype
px.11  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.11(Fii, Gii, par.list)* Wm.mat[1,1]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.13  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.13(Fii, Gii, par.list)* Wm.mat[1,3]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.33  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.33(Fii, Gii, par.list)* Wm.mat[3,3]) / PxTot(Fii, Gii, par.list, Wm.mat)
}

# Linkage Disequilibrium (convenience function for outcross haplotype frequency equations)
LD  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	LD  <-  par.list$r * ((ox.14(Fii, Gii, par.list, Wf.mat) - ox.23(Fii, Gii, par.list, Wf.mat)) / 2)
	if(is.nan(LD)) 0 
		else(LD)
}

#########################################
## Haplotype frequency change among outcrossed gametes

# Ovules
x1  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(ox.11(Fii, Gii, par.list, Wf.mat) + ((ox.12(Fii, Gii, par.list, Wf.mat) + ox.13(Fii, Gii, par.list, Wf.mat) + ox.14(Fii, Gii, par.list, Wf.mat))/2)) - LD(Fii, Gii, par.list, Wf.mat)
}
x2  <-  function(Fii, Gii, par.list, Wf.mat, ...){
	(ox.22(Fii, Gii, par.list, Wf.mat) + ((ox.12(Fii, Gii, par.list, Wf.mat) + ox.23(Fii, Gii, par.list, Wf.mat) + ox.24(Fii, Gii, par.list, Wf.mat))/2)) + LD(Fii, Gii, par.list, Wf.mat)
}
x3  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(ox.33(Fii, Gii, par.list, Wf.mat) + ((ox.13(Fii, Gii, par.list, Wf.mat) + ox.23(Fii, Gii, par.list, Wf.mat) + ox.34(Fii, Gii, par.list, Wf.mat))/2)) + LD(Fii, Gii, par.list, Wf.mat)
}
x4  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(ox.44(Fii, Gii, par.list, Wf.mat) + ((ox.14(Fii, Gii, par.list, Wf.mat) + ox.24(Fii, Gii, par.list, Wf.mat) + ox.34(Fii, Gii, par.list, Wf.mat))/2)) - LD(Fii, Gii, par.list, Wf.mat)
}

# Pollen/Sperm
y1  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	px.11(Fii, Gii, par.list, Wm.mat) + (px.13(Fii, Gii, par.list, Wm.mat) / 2)
}
y2  <-  function(Fii, Gii, par.list, Wm.mat, ...){
	0
}
y3  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	 	px.33(Fii, Gii, par.list, Wm.mat) + (px.13(Fii, Gii, par.list, Wm.mat) / 2)
}
y4  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	0
}


#######################################################
## Genotypic frequency x transmission mode recursions
#######################################################

# Genotypic frequency recursions for offspring produced by outcross fertilization
Fpr.11  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x1(Fii, Gii, par.list, Wf.mat)*y1(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.12  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x2(Fii, Gii, par.list, Wf.mat)*y1(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.13  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x1(Fii, Gii, par.list, Wf.mat)*y3(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.14  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x4(Fii, Gii, par.list, Wf.mat)*y1(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.22  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	0
}
Fpr.23  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x2(Fii, Gii, par.list, Wf.mat)*y3(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.24  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	0
}
Fpr.33  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x3(Fii, Gii, par.list, Wf.mat)*y3(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.34  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x4(Fii, Gii, par.list, Wf.mat)*y3(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.44  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
    0
}

# Genotypic frequency recursions for offspring produced by self-fertlization
Gpr.11  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(os.11(Fii, Gii, par.list, Wf.mat) + (os.13(Fii, Gii, par.list, Wf.mat)/4)) * S(Fii, Gii, par.list, Wf.mat);
}
Gpr.12  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	0
}
Gpr.13  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(os.13(Fii, Gii, par.list, Wf.mat)/2) * S(Fii, Gii, par.list, Wf.mat);
}
Gpr.14  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	0
}
Gpr.22  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	0
}
Gpr.23  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	0
}
Gpr.24  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	0
}
Gpr.33  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(os.33(Fii, Gii, par.list, Wf.mat) + (os.13(Fii, Gii, par.list, Wf.mat)/4)) * S(Fii, Gii, par.list, Wf.mat);
}
Gpr.34  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	0
}
Gpr.44  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
    0
}


# Single-locus sterility allele equilibrium frequencies, Eq.5 in Charlesworth (1978).
Zhat.gyn  <-  function(par.list) {
	(par.list$k + 2*par.list$C*par.list$delta - 1) / (2*(par.list$k + par.list$C*par.list$delta))
}



################################################################################################
################################################################################################


#' Forward deterministic simulation of genotypic recursions for
#' 2-locus SA w/ linked sterility locus
#'
#' @title Forward deterministic simulation of genotypic recursions
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'				   gen    =  5000,
#'				   C      =  0,
#'                 delta  =  0,
#'				   sm     =  0.7,
#'				   sf     =  0.7,
#'				   hm     =  1/2,
#'				   hf     =  1/2,
#'                 k      =  1,
#'				   r      =  0.5
#'				   )
#' @param Fii.init A vector of initial genotypic frequencies for offspring produced via outcrossing (must have length = 10).
#' c((1 - C) freqAA,0,(1 - C) freqAa,0,0,0,0,(1 - C) freqaa,0,0) for invasion of aaM1M2 into population 'fixed' for AABB.
#' c((1 - C) freqAA,0,(1 - C) freqAa,0,0,0,0,(1 - C) freqaa,0,0) for invasion of AAM1M2 into population 'fixed' for aabb.
#' @param Gii.init A vector of initial genotypic frequencies for offspring produced via self-fertilization (must have length = 10).
#' c(C*freqAA,0,C*freqAa,0,0,0,0,C*freqaa,0,0) for invasion of aaM1M2 into population 'fixed' for AABB.
#' c(C*freqAA,0,C*freqAa,0,0,0,0,C*freqaa,0,0) for invasion of AAM1M2 into population 'fixed' for aabb.
#' @return Returns a list with timeseries for each genotype, equilibrium frequencies, and a numeric (0,1) for whether the 
#' equilibrium was polymorphic (with tolerance 1E-6).
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSim(par.list, Fii.init, Gii.init, threshold = 1e-6) 
recursionFwdSim  <-  function(par.list, Fii.init, Gii.init, threshold = 1e-6, ...) {

	##  Warnings
	if(any(par.list[2:7] < 0) | any(par.list[2:7] > 1) | par.list$r > 0.5)
		stop('The chosen parameter values fall outside of the reasonable bounds')

	if(par.list$hf  !=  par.list$hm)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(par.list$hf != 0.5 & par.list$hf != 0.25)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(sum(Fii.init, Gii.init) != 1)
		stop('Incorrect initial frequencies. Initial frequencies must sum to 1')

	##  Fitness Expression Matrices
	Wf.mat   <-  matrix(
						c(1, (1 + par.list$k), (1 - par.list$hf*par.list$sf),                  (1 - par.list$hf*par.list$sf)*(1 + par.list$k),
                          0, (1 + par.list$k), (1 - par.list$hf*par.list$sf)*(1 + par.list$k), (1 - par.list$hf*par.list$sf)*(1 + par.list$k),
                          0, 0,                (1 - par.list$sf),                              (1 - par.list$sf)*(1 + par.list$k),
                          0, 0,                0,                                              (1 - par.list$sf)*(1 + par.list$k)), 
						nrow=4, byrow=TRUE
						)

	Wm.mat   <-  matrix(
						c((1 - par.list$sm), 0, (1 - par.list$hm*par.list$sm), 0,
                          0,                 0, 0,                             0,
                          0,                 0, (1 - par.list$sm),             0,
                          0,                 0, 0,                             0), 
						nrow=4, byrow=TRUE
						)	


	##  Initilize data storage structures
	Fii.gen  <-  matrix(0, ncol=20, nrow=par.list$gen)
	colnames(Fii.gen)  <-  c('Fpr.11', 'Fpr.12', 'Fpr.13', 'Fpr.14', 'Fpr.22', 'Fpr.23', 'Fpr.24', 'Fpr.33', 'Fpr.34', 'Fpr.44',
		                     'Gpr.11', 'Gpr.12', 'Gpr.13', 'Gpr.14', 'Gpr.22', 'Gpr.23', 'Gpr.24', 'Gpr.33', 'Gpr.34', 'Gpr.44')
	
	##  Generation Loop
		# initialize
	for (j in 1:ncol(Fii.gen)) {
		recFct          <-  get(colnames(Fii.gen)[j])
		Fii.gen[1, j]   <-  round(recFct(Fii = Fii.init, Gii = Gii.init, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat), digits=6)
	}

	# Start simulation
	i      <-  2
	diffs  <-  rep(1,20)
	while (i < par.list$gen & any(diffs[diffs != 0] > threshold)) {
		for (j in 1:ncol(Fii.gen)) {
			recFct          <-  get(colnames(Fii.gen)[j])
			Fii.gen[i, j]   <-  round(recFct(Fii = Fii.gen[i-1, 1:10], Gii = Fii.gen[i-1, 11:20], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		}
		
		diffs  <-  Fii.gen[i,] - Fii.gen[i-1,]
		i      <-  i+1
	}

	# Calculate 1-locus equilibrium frequency of unisexuals
	Zhat  <-  Zhat.gyn(par.list)
	qHat  <-  qHatAdd(C = par.list$C, delta = par.list$delta, sf = par.list$sf, sm = par.list$sm)

	##  Output
	res  <-  list(
				  "par.list" =  par.list,
				  "Fii.gen"  =  Fii.gen[1:i-1,],
				  "EQ.freq"  =  Fii.gen[i-1,],
				  "Zhat"     =  Zhat,
				  "qHat"     =  qHat
 				 )
	return(res)
}




#' Simulation loop wrapping forward deterministic simulations 
#' of genotypic recursions 
#' #'
#' @title Forward deterministic simulation of genotypic recursions.
#' @param n number of randomly generated values for sf & sm. 
#' Determines resolution with which parameter space is explored.
#' @param gen Maximum number of generations for each simulation (as in par.list).
#' @param C The fixed selfing rate (as in par.list).
#' @param hf Dominance through female expression (as in par.list).
#' @param hm Dominance through male expression (as in par.list).
#' @param r.vals Values of recombination rate to explore(as in par.list).
#' @param threshold Threshold difference between genotypic frequencies before simulation cuts off.
#' @return Returns a data frame with parameter values, a variable describing whether 
#' the final state of the simulation was polymorphic polymorphism, whether evaluating the eigenvalues
#' predicts polymorphism, and whether these two methods agree with one another.
#' @seealso `recursionFwdSim`
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0, delta = 0, 
#'                     hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
#'                     seed = 3497016, threshold = 1e-7)
recursionFwdSimLoop  <-  function(gen = 5000, C = 0, delta = 0, hf = 0.5, hm = 0.5, 
	                              sm.vals = c(0.05, 0.2, 0.8, 0.95), r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
	                              threshold = 1e-7) {

	## Warnings
	if(any(c(C,hf,hm,r.vals) < 0) | any(c(C,hf,hm) > 1) | any(r.vals > 0.5))
		stop('At least one of the chosen parameter values fall outside of the reasonable bounds')

	if(threshold > 1e-7)
		stop('Carefully consider whether you want to change this threshold, 
			  as it will affect how many generations are required to reach
			  convergence, and thus how long the simulations take')

	# Calculate sf values
	s.vals  <-  c(0,0)
	for (i in 1:length(sm.vals)) {
			s.vals  <-  rbind(s.vals, cbind(funnelSlice(sm.vals[i], C, delta, h = hf), sm.vals[i]))
		}
	s.vals  <-  s.vals[-1,]

	# Calculate k values
	Ks  <-  invGyn(C, delta) + c(0.2, 0.1, 0, -0.1, -0.2)

	#  initialize selection coeficients and storage structures
	eqFreqs  <-  matrix(0, nrow=length(r.vals)*length(Ks)*nrow(s.vals), ncol=20)
	qHat     <-  rep(0, length(r.vals)*length(Ks)*nrow(s.vals))
	
	##  Simulation Loop over values of r, sm, sf for fixed selfing rate (C)
	for (i in 1:length(r.vals)) {
		for (j in 1:length(Ks)) {
			for (m in 1:nrow(s.vals)) {
				
				par.list  <-  list(
							   		gen    =  5000,
									C      =  C,
									delta  =  delta,
									sf     =  s.vals[m,1],
									sm     =  s.vals[m,2],
									hf     =  hm,
									hm     =  hf,
				            		k      =  Ks[j],
									r      =  r.vals[i]
								   )

				##  Set initial frequencies to speed up convergence. 
				if(par.list$sf > Inv.a.add(par.list$sm, par.list$C, par.list$delta)) {
					Fii.init    <-  c((1 - C)*0.98,0,(1 - C)*0.01,(1 - C)*0.01,0,0,0,(1 - C)*0,0,0)
					Gii.init    <-  c((C*0.98),0,(C*0.01),(C*0.01),0,0,0,(C*0),0,0)
				}
				if(par.list$sf < Inv.A.add(par.list$sm, par.list$C, par.list$delta)) {
					Fii.init    <-  c((1 - C)*0,0,(1 - C)*0.01,(1 - C)*0.01,0,0,0,(1 - C)*0.98,0,0)
					Gii.init    <-  c((C*0),0,(C*0.01),(C*0.01),0,0,0,(C*0.98),0,0)
				}
				if(par.list$sf < Inv.a.add(par.list$sm, par.list$C, par.list$delta) & 
					   par.list$sf > Inv.A.add(par.list$sm, par.list$C, par.list$delta)) {
					if(hm == 0.5 & hf == 0.5) {
					   	qhat  <-  qHatAdd(par.list$C, par.list$delta, par.list$sf, par.list$sm)
					   	QEs   <-  c(QE.FAA(q = qhat, C = par.list$C), QE.FAa(q = qhat, C = par.list$C), QE.Faa(q = qhat, C = par.list$C))
					   	QEs  <- QEs - 0.01
						Fii.init    <-  as.numeric(rounded(c((1 - C)*QEs[1], (1 - C)*0.01, (1 - C)*QEs[2], (1 - C)*0.01, 0, 0, 0, (1 - C)*QEs[3], (1 - C)*0.01, 0), precision = 15))
						Gii.init    <-  as.numeric(rounded(c(      C*QEs[1],       C*0.01,       C*QEs[2],       C*0.01, 0, 0, 0,       C*QEs[3],       C*0.01, 0), precision = 15))
					}
					if(hf == hm & hm < 0.5 & hf < 0.5) {
					   	qhat  <-  qHatDomRev(par.list$C, par.list$delta, par.list$sf, par.list$sm, par.list$hf)
					   	QEs   <-  c(QE.FAA(q = qhat, C = par.list$C), QE.FAa(q = qhat, C = par.list$C), QE.Faa(q = qhat, C = par.list$C))
					   	QEs  <- QEs - 0.01
						Fii.init    <-  as.numeric(rounded(c((1 - C)*QEs[1], (1 - C)*0.01, (1 - C)*QEs[2], (1 - C)*0.01, 0, 0, 0, (1 - C)*QEs[3], (1 - C)*0.01, 0), precision = 15))
						Gii.init    <-  as.numeric(rounded(c(      C*QEs[1],       C*0.01,       C*QEs[2],       C*0.01, 0, 0, 0,       C*QEs[3],       C*0.01, 0), precision = 15))
					}
				}

				# Run simulation for given parameter values
				res      <-  recursionFwdSim(par.list = par.list, Fii.init = Fii.init, Gii.init = Gii.init, threshold = threshold)

				# Store equilibrium frequencies
				eqFreqs[((i-1)*length(Ks)*nrow(s.vals)) + ((j-1)*nrow(s.vals)) + m,]  <-  res$EQ.freq
				qHat[((i-1)*length(Ks)*nrow(s.vals))    + ((j-1)*nrow(s.vals)) + m]   <-  res$qHat
			}
		}
	}

	#  Compile results as data.frame
	rs  <-  c()
	ks  <-  c()
	for(i in 1:length(r.vals)) {
		rs  <-  c(rs, rep(r.vals[i], length(Ks)*nrow(s.vals)))
	}
	for(i in 1:length(r.vals)) {
		for(j in 1:length(Ks)) {
			ks  <-  c(ks, rep(Ks[j], nrow(s.vals)))
		}
	}
	results.df  <-  data.frame("hf"      =  rep(hf, length(r.vals)*length(Ks)*nrow(s.vals)),
							   "hm"      =  rep(hm, length(r.vals)*length(Ks)*nrow(s.vals)),
							   "C"       =  rep(C,   length(r.vals)*length(Ks)*nrow(s.vals)),
							   "delta"   =  rep(delta,   length(r.vals)*length(Ks)*nrow(s.vals)),
							   "r"       =  rs,
							   "k"       =  ks,
							   "sf"      =  rep(s.vals[,1], length(r.vals)*length(Ks)),
							   "sm"      =  rep(s.vals[,2], length(r.vals)*length(Ks))
							   )
	results.df  <-  cbind(results.df, eqFreqs, qHat)
	colnames(results.df)  <-  c('hf','hm','C','delta','r','k','sf','sm',
		                        'F.11', 'F.12', 'F.13', 'F.14', 'F.22', 'F.23', 'F.24', 'F.33', 'F.34', 'F.44', 
	    	                    'G.11', 'G.12', 'G.13', 'G.14', 'G.22', 'G.23', 'G.24', 'G.33', 'G.34', 'G.44', 'qHat')
	
	#  Write results.df to .txt file
	filename  <-  paste("./output/data/simResults/gyn-dom", "_C", C, "_delta", delta, "_add", "_strgSel", ".csv", sep="")
	write.csv(results.df, file=filename, row.names = FALSE)

	#  Return results.df in case user wants it
	return(results.df)
}
