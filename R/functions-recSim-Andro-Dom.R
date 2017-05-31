#####################################################
#  Androdioecy via invasion of dominant 
#  female-sterility mutations
#
#  Necessary functions for deterministic simulations
#  of genotypic frequency recursions for a model of
#  the evolution of andodioecy via invasion of a
#  completely dominant female sterility allele 
#
#  Author: Colin Olito
#
#  NOTES:  
#          

##########################
##  Key to Fii subscripts
#  F11  =  Fii[1]   =  AAM1M1
#  F12  =  Fii[2]   =  AAM1M2
#  F13  =  Fii[3]   =  AaM1M1
#  F14  =  Fii[4]   =  AaM1M2
#  F22  =  Fii[5]   =  AAM2M2
#  F23  =  Fii[6]   =  aAM2M1
#  F24  =  Fii[7]   =  AaM2M2
#  F33  =  Fii[8]   =  aaM1M1
#  F34  =  Fii[9]   =  aaM1M2
#  F44  =  Fii[10]  =  aaM2M2

#######################################################
## Necessary Functions
#######################################################

# Average fitness of adults after inbreeding depression, but prior to mating
Dbar  <-  function(Fii, Gii, par.list, ...) {
	(Fii[1] + Fii[2] + Fii[3] + Fii[4] + Fii[5] + Fii[6] + Fii[7] + Fii[8] + Fii[9] + Fii[10]) + (1 - par.list$delta)*(Gii[1] + Gii[2] + Gii[3] + Gii[4] + Gii[5] + Gii[6] + Gii[7] + Gii[8] + Gii[9] + Gii[10])
}

# Frequencies of adult genotypes after inbreeding depression 
FA.11  <-  function(Fii, Gii, par.list, ...) {
	(Fii[1] + Gii[1]   * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.12  <-  function(Fii, Gii, par.list, ...) {
	(Fii[2] + Gii[2]   * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.13  <-  function(Fii, Gii, par.list, ...) {
	(Fii[3] + Gii[3]   * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.14  <-  function(Fii, Gii, par.list, ...) {
	(Fii[4] + Gii[4]   * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.22  <-  function(Fii, Gii, par.list, ...) {
	(Fii[5] + Gii[5]   * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.23  <-  function(Fii, Gii, par.list, ...) {
	(Fii[6] + Gii[6]   * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.24  <-  function(Fii, Gii, par.list, ...) {
	(Fii[7] + Gii[7]   * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.33  <-  function(Fii, Gii, par.list, ...) {
	(Fii[8] + Gii[8]   * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.34  <-  function(Fii, Gii, par.list, ...) {
	(Fii[9] + Gii[9]   * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}
FA.44  <-  function(Fii, Gii, par.list, ...) {
	(Fii[10] + Gii[10] * (1 - par.list$delta)) / Dbar(Fii, Gii, par.list)
}

# Total ovules produced 
OTot  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	FA.11(Fii, Gii, par.list)*Wf.mat[1,1] + FA.12(Fii, Gii, par.list)*Wf.mat[1,2] + FA.13(Fii, Gii, par.list)*Wf.mat[1,3]
}

# Total ovule used for self-fertilization
OsTot  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	par.list$C * OTot(Fii, Gii, par.list, Wf.mat)
}

# Total ovules used for outcrossing
OxTot  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
		(1 - par.list$C)*(FA.11(Fii, Gii, par.list)*Wf.mat[1,1] + FA.13(Fii, Gii, par.list)*Wf.mat[1,3] + FA.33(Fii, Gii, par.list)*Wf.mat[3,3])
}

# Proportion of all ovules that are self-fertilized
S  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	OsTot(Fii, Gii, par.list, Wf.mat) / OTot(Fii, Gii, par.list, Wf.mat)
}

# Proportional contribution to selfed offspring for each genotype capable of selfing
os.11  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	if(OsTot(Fii, Gii, par.list, Wf.mat) == 0) 0
		else(par.list$C*(FA.11(Fii, Gii, par.list)*Wf.mat[1,1]) / OsTot(Fii, Gii, par.list, Wf.mat))
}
os.13  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	if(OsTot(Fii, Gii, par.list, Wf.mat) == 0) 0
		else(par.list$C*(FA.13(Fii, Gii, par.list)*Wf.mat[1,3]) / OsTot(Fii, Gii, par.list, Wf.mat))
}
os.33  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	if(OsTot(Fii, Gii, par.list, Wf.mat) == 0) 0
		else(par.list$C*(FA.33(Fii, Gii, par.list)*Wf.mat[3,3]) / OsTot(Fii, Gii, par.list, Wf.mat))
}

# Proportional contribution to outcrossed offspring for each genotype capable of outcrossing
ox.11  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(1 - par.list$C)*(FA.11(Fii, Gii, par.list)*Wf.mat[1,1]) / OxTot(Fii, Gii, par.list, Wf.mat)
}
ox.13  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(1 - par.list$C)*(FA.13(Fii, Gii, par.list)*Wf.mat[1,3]) / OxTot(Fii, Gii, par.list, Wf.mat)
}
ox.33  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	(1 - par.list$C)*(FA.33(Fii, Gii, par.list)*Wf.mat[3,3]) / OxTot(Fii, Gii, par.list, Wf.mat)
}

# Total pollen used for outcrossing (proportional to) the total amont of pollen 
# in the population (the amount of pollen used for selfing has negligable effect 
# on the pool of outcrossing pollen)
PxTot  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.11(Fii, Gii, par.list)* Wm.mat[1,1] + FA.13(Fii, Gii, par.list)* Wm.mat[1,3] + FA.33(Fii, Gii, par.list)* Wm.mat[3,3]) + 
	(FA.12(Fii, Gii, par.list)*Wm.mat[1,2] + FA.14(Fii, Gii, par.list)*Wm.mat[1,4] + FA.22(Fii, Gii, par.list)*Wm.mat[2,2] + 
	    FA.23(Fii, Gii, par.list)*Wm.mat[2,3] + FA.24(Fii, Gii, par.list)*Wm.mat[2,4] + FA.34(Fii, Gii, par.list)*Wm.mat[3,4] + 
        FA.44(Fii, Gii, par.list)*Wm.mat[4,4])
}

# Frequency of pollen/sperm | Outcrossing genotype
px.11  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.11(Fii, Gii, par.list)* Wm.mat[1,1]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.12  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.12(Fii, Gii, par.list)* Wm.mat[1,2]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.13  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.13(Fii, Gii, par.list)* Wm.mat[1,3]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.14  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.14(Fii, Gii, par.list)* Wm.mat[1,4]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.22  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.22(Fii, Gii, par.list)* Wm.mat[2,2]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.23  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.23(Fii, Gii, par.list)* Wm.mat[2,3]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.24  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.24(Fii, Gii, par.list)* Wm.mat[2,4]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.33  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.33(Fii, Gii, par.list)* Wm.mat[3,3]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.34  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.34(Fii, Gii, par.list)* Wm.mat[3,4]) / PxTot(Fii, Gii, par.list, Wm.mat)
}
px.44  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	(FA.44(Fii, Gii, par.list)* Wm.mat[4,4]) / PxTot(Fii, Gii, par.list, Wm.mat)
}

# Linkage Disequilibrium (convenience function for outcross haplotype frequency equations)
LD  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	par.list$r*((px.14(Fii, Gii, par.list, Wm.mat) - px.23(Fii, Gii, par.list, Wm.mat)) / 2)
}

#########################################
## Haplotype frequency change among outcrossed gametes

# Ovules
x1  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	ox.11(Fii, Gii, par.list, Wf.mat) + (ox.13(Fii, Gii, par.list, Wf.mat) / 2)
}
x2  <-  function(Fii, Gii, par.list, Wf.mat, ...){
	0
}
x3  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	 	ox.33(Fii, Gii, par.list, Wf.mat) + (ox.13(Fii, Gii, par.list, Wf.mat) / 2)
}
x4  <-  function(Fii, Gii, par.list, Wf.mat, ...) {
	0
}

# Pollen/Sperm
y1  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	((2*px.11(Fii, Gii, par.list, Wm.mat) + px.12(Fii, Gii, par.list, Wm.mat) + px.13(Fii, Gii, par.list, Wm.mat) + px.14(Fii, Gii, par.list, Wm.mat))/2) - LD(Fii, Gii, par.list, Wm.mat)
}
y2  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	((2*px.22(Fii, Gii, par.list, Wm.mat) + px.12(Fii, Gii, par.list, Wm.mat) + px.23(Fii, Gii, par.list, Wm.mat) + px.24(Fii, Gii, par.list, Wm.mat))/2) + LD(Fii, Gii, par.list, Wm.mat)
}
y3  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	((2*px.33(Fii, Gii, par.list, Wm.mat) + px.13(Fii, Gii, par.list, Wm.mat) + px.23(Fii, Gii, par.list, Wm.mat) + px.34(Fii, Gii, par.list, Wm.mat))/2) + LD(Fii, Gii, par.list, Wm.mat)
}
y4  <-  function(Fii, Gii, par.list, Wm.mat, ...) {
	((2*px.44(Fii, Gii, par.list, Wm.mat) + px.14(Fii, Gii, par.list, Wm.mat) + px.24(Fii, Gii, par.list, Wm.mat) + px.34(Fii, Gii, par.list, Wm.mat))/2) - LD(Fii, Gii, par.list, Wm.mat)
}


#######################################################
## Genotypic frequency x transmission mode recursions
#######################################################

# Genotypic frequency recursions for offspring produced by outcross fertilization
Fpr.11  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x1(Fii, Gii, par.list, Wf.mat)*y1(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.12  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x1(Fii, Gii, par.list, Wf.mat)*y2(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.13  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x1(Fii, Gii, par.list, Wf.mat)*y3(Fii, Gii, par.list, Wm.mat) + 
     x3(Fii, Gii, par.list, Wf.mat)*y1(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.14  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x1(Fii, Gii, par.list, Wf.mat)*y4(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.22  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	0
}
Fpr.23  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x3(Fii, Gii, par.list, Wf.mat)*y2(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.24  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	0
}
Fpr.33  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x3(Fii, Gii, par.list, Wf.mat)*y3(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
}
Fpr.34  <-  function(Fii, Gii, par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat,...) {
	(x3(Fii, Gii, par.list, Wf.mat)*y4(Fii, Gii, par.list, Wm.mat))*(1 - S(Fii, Gii, par.list, Wf.mat));
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
recursionFwdSim  <-  function(par.list, Fii.init, Gii.init, threshold = 1e-7, ...) {

	##  Warnings
	if(any(par.list[2:7] < 0) | any(par.list[2:7] > 1) | par.list$r > 0.5)
		stop('The chosen parameter values fall outside of the reasonable bounds')

	if(par.list$hf  !=  par.list$hm)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(par.list$hf != 0.5 & par.list$hf != 0.25)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(round(sum(Fii.init, Gii.init), digits=4) != 1)
		stop('Incorrect initial frequencies. Initial frequencies must sum to 1')

	##  Fitness Expression Matrices
	Wf.mat   <-  matrix(
						c(1, 0, (1 - par.list$hf*par.list$sf), 0,
                          0, 0,                             0, 0,
                          0, 0,             (1 - par.list$sf), 0,
                          0, 0,                             0, 0), 
                         nrow=4, byrow=TRUE
						)

	Wm.mat   <-  matrix(
						c((1 - par.list$sm), (1 - par.list$sm)*(1 + par.list$k), (1 - par.list$hm*par.list$sm),                  (1 - par.list$hm*par.list$sm)*(1 + par.list$k),
                                          0, (1 - par.list$sm)*(1 + par.list$k), (1 - par.list$hm*par.list$sm)*(1 + par.list$k), (1 - par.list$hm*par.list$sm)*(1 + par.list$k),
                                          0,                                  0,                                              1,                               (1 + par.list$k),
                                          0,                                  0,                                              0,                               (1 + par.list$k)), 
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
			Fii.gen[i, j]   <-  round(recFct(Fii = Fii.gen[i-1, 1:10], Gii = Fii.gen[i-1, 11:20], par.list = par.list, Wf.mat = Wf.mat, Wm.mat = Wm.mat), digits=6)
		}

		diffs  <-  abs(Fii.gen[i,] - Fii.gen[i-1,])
		i      <-  i+1
	}
	if(i == par.list$gen) {
		print('Warning: maximum runtime reached. Results may not represent equilibrium frequencies')
	}

	# Calculate 1-locus equilibrium frequency of unisexuals
	Zhat  <-  Zhat.and.list(par.list)
	if(par.list$hf == par.list$hf & par.list$hf == 0.5) {
		qHat  <-  qHatAdd(C = par.list$C, delta = par.list$delta, sf = par.list$sf, sm = par.list$sm)
	}
	if(par.list$hf == par.list$hf & par.list$hf == 0.25) {
		qHat  <-  qHatDomRev(C = par.list$C, delta = par.list$delta, sf = par.list$sf, sm = par.list$sm)
	}
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
#' 
#' @title Forward deterministic simulation of genotypic recursions.
#' @param gen Maximum number of generations for each simulation (as in par.list).
#' @param dStar Determines the maximum level of inbreeding depression (see `deltaC` function)
#' @param a Determines linearity of `deltaC` results
#' @param b Determines range of inbreeding depression from `deltaC` results
#' @param sf Selection through female expression (as in par.list).
#' @param sm Selection through male expression (as in par.list).
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
recursionFwdSimLoop  <-  function(gen = 5000, dStar = 0.8, a = 1, b = 0.5, 
	                              sm = 0.5, hf = 0.5, hm = 0.5, 
	                              kMult = c(1.1, 0.95, 0.9), r.vals = c(0.0, 0.01, 0.05, 0.1), 
	                              threshold = 1e-7) {

	## Warnings
	if(any(c(dStar,sm,hf,hm,r.vals) < 0) | any(c(dStar,sm,hf,hm) > 1) | any(r.vals > 0.5))
		stop('At least one of the chosen parameter values fall outside of the reasonable bounds')

	if(threshold > 1e-6)
		stop('Carefully consider whether you want to change this threshold, 
			  as it will affect how many generations are required to reach
			  convergence, and thus how long the simulations take')

	# Calculate Cs, deltas, sfs 
	Cs  <-  seq(0,0.9,by=0.1)
	Ds  <-  deltaC(dStar=dStar, C=Cs, a=a, b=b)

	if(hf == hm & hm == 0.5) {
		sfs  <-  equalPQ.add(sm = sm, C = Cs, delta = Ds)
	}
	if(hf == hm & hm == 0.25) {
		sfs  <-  equalPQ.domRev(sm = sm, C = Cs, delta = Ds)	
	}

	# Calculate k values
	Ks     <-  matrix(0, nrow=length(kMult), ncol=length(Cs))
	for(i in 1:nrow(Ks)){
    	Ks[i,]  <-  invAnd(Cs, Ds) * kMult[i]
    }

	#  initialize storage structures
	eqFreqs  <-  matrix(0, nrow=length(r.vals)*nrow(Ks)*length(Cs), ncol=20)
	qHat     <-  rep(0, length(r.vals)*nrow(Ks)*length(Cs))
	ZHat     <-  rep(0, length(r.vals)*nrow(Ks)*length(Cs))
	
	##  Simulation Loop over values of r, sm, sf for fixed selfing rate (C)
	print('Running Deterministic Recursion Simulations')
	pb   <-  txtProgressBar(min=0, max=length(r.vals)*nrow(Ks), style=3)
	setTxtProgressBar(pb, 0)
	for (i in 1:length(r.vals)) {
		for (j in 1:nrow(Ks)) {
			for (m in 1:length(Cs)) {
				
				par.list  <-  list(
							   		gen    =  gen,
									C      =  Cs[m],
									delta  =  Ds[m],
									sf     =  sfs[m],
									sm     =  sm,
									hf     =  hm,
									hm     =  hf,
									k      =  Ks[j,m],
									r      =  r.vals[i]
								   )

				##  Set initial frequencies 
					# Additive effects
					if(hm == 0.5 & hf == 0.5) {
					   	qhat  <-  qHatAdd(par.list$C, par.list$delta, par.list$sf, par.list$sm)
					   }
					# Dominance Reversal
					if(hm == 0.25 & hf == 0.25) {
					   	qhat  <-  qHatDomRev(par.list$C, par.list$delta, par.list$sf, par.list$sm, par.list$hf)
					   	}
				   	QEs  <-  c(QE.FAA(q = qhat, C = par.list$C), QE.FAa(q = qhat, C = par.list$C), QE.Faa(q = qhat, C = par.list$C))
				   	if(par.list$C == 0) {
						QEs[QEs == max(QEs)]  <-  max(QEs) - 0.01
						Fii.init  <-  round(c(QEs[1], 0, QEs[2], 0.01, 0, 0, 0, QEs[3], 0, 0), digits=5)
						Gii.init  <-  rep(0,10)
				   	}
				   	else {
				   		QEs   <-  QEs/sum(QEs)
				   	   	QEs[QEs == max(QEs)][1]  <-  max(QEs) - 0.01
						Fii.init  <-  round(c((1 - par.list$C)*QEs[1], 0, (1 - par.list$C)*QEs[2], (1 - par.list$C)*0.01, 0, 0, 0, (1 - par.list$C)*QEs[3], 0, 0), digits=5)
						Gii.init  <-  round(c(par.list$C*QEs[1], 0, par.list$C*QEs[2], par.list$C*0.01, 0, 0, 0, par.list$C*QEs[3], 0, 0), digits=5)
					}
					if(sum(Fii.init,Gii.init) != 1) {
						Fii.init  <-  Fii.init/sum(Fii.init, Gii.init)
						Gii.init  <-  Gii.init/sum(Fii.init, Gii.init)
					}

		 		# Run simulation for given parameter values
				res  <-  recursionFwdSim(par.list = par.list, Fii.init = Fii.init, Gii.init = Gii.init, threshold = threshold)

				# Store equilibrium frequencies
				eqFreqs[((i-1)*nrow(Ks)*length(Cs)) + ((j-1)*length(Cs)) + m,]  <-  res$EQ.freq
				qHat[((i-1)*nrow(Ks)*length(Cs))    + ((j-1)*length(Cs)) + m]   <-  res$qHat
				ZHat[((i-1)*nrow(Ks)*length(Cs))    + ((j-1)*length(Cs)) + m]   <-  res$Zhat
			}
		setTxtProgressBar(pb, ((i-1)*nrow(Ks) + j))
		}
	}
	setTxtProgressBar(pb, length(r.vals)*nrow(Ks))

	#  Compile results as data.frame
	rs  <-  c()
	ks  <-  c()
	for(i in 1:length(r.vals)) {
		rs  <-  c(rs, rep(r.vals[i], nrow(Ks)*length(Cs)))
	}
	for(i in 1:length(r.vals)) {
		for(j in 1:nrow(Ks)) {
			ks  <-  c(ks, rep(kMult[j], length(Cs)))
		}
	}
	results.df  <-  data.frame("hf"      =  rep(hf, length(r.vals)*nrow(Ks)*length(Cs)),
							   "hm"      =  rep(hm, length(r.vals)*nrow(Ks)*length(Cs)),
							   "C"       =  rep(Cs, length(r.vals)*nrow(Ks)),
							   "delta"   =  rep(Ds, length(r.vals)*nrow(Ks)),
							   "r"       =  rs,
							   "k"       =  ks,
							   "sf"      =  rep(sfs, length(r.vals)*nrow(Ks)),
							   "sm"      =  sm
							   )
	results.df  <-  cbind(results.df, eqFreqs, qHat,ZHat)
	colnames(results.df)  <-  c('hf','hm','C','delta','r','k','sf','sm',
		                        'F.11', 'F.12', 'F.13', 'F.14', 'F.22', 'F.23', 'F.24', 'F.33', 'F.34', 'F.44', 
	    	                    'G.11', 'G.12', 'G.13', 'G.14', 'G.22', 'G.23', 'G.24', 'G.33', 'G.34', 'G.44', 'qHat','ZHat')
	
	#  Write results.df to .txt file
	if(hm == hf & hf == 0.5) {
		filename  <-  paste("./output/data/simResults/and-dom", "_dStar", dStar, "_a", a, "_sm", sm, "_add", ".csv", sep="")
	}
	if(hm == hf & hf == 0.25) {
		filename  <-  paste("./output/data/simResults/and-dom", "_dStar", dStar, "_a", a, "_sm", sm, "_domRev", ".csv", sep="")
	}
	write.csv(results.df, file=filename, row.names = FALSE)

	#  Return results.df in case user wants it
	return(results.df)
}





#######################################
##  Plotting functions

#' Exploratory plots: Equilibrium frequencies of 
#'                    -- Additive fitness effects at SA locus
#' @title Invasion of dominant male sterility allele into populations
#' @param df .csv file with data output from recursion simulations.
#' @examples andDomRecPlots(df = "./output/data/simResults/and-dom_dStar0.8_a1_sm0.5_add.csv")
#' @author Colin Olito
#' @export
#' 
andDomRecPlots  <-  function(df = "./output/data/simResults/and-dom_dStar0.8_a1_sm0.4_add.csv") {

    # Import data
    data  <-  read.csv(df, header=TRUE)

	# Calculate equilibrium frequencies of M2, males, A, a
	data$q.m2        <-  (data$'F.12'/2) + (data$'F.14'/2) + data$'F.22' + (data$'F.23'/2) + data$'F.24' + (data$'F.34'/2) + data$'F.44' + 
	                     (data$'G.12'/2) + (data$'G.14'/2) + data$'G.22' + (data$'G.23'/2) + data$'G.24' + (data$'G.34'/2) + data$'G.44'
	data$males     <-  data$'F.12' + data$'F.14' + data$'F.22' + data$'F.23' + data$'F.24' + data$'F.34' + data$'F.44' + 
	                     data$'G.12' + data$'G.14' + data$'G.22' + data$'G.23' + data$'G.24' + data$'G.34' + data$'G.44'
	data$p.A         <-  data$'F.11' + data$'F.12' + (data$'F.13'/2) + (data$'F.14'/2) + data$'F.22' + (data$'F.23'/2) + (data$'F.24'/2) + 
	                     data$'G.11' + data$'G.12' + (data$'G.13'/2) + (data$'G.14'/2) + data$'G.22' + (data$'G.23'/2) + (data$'G.24'/2)
	data$q.a         <-  (data$'F.13'/2) + (data$'F.14'/2) + (data$'F.23'/2) + (data$'F.24'/2) + data$'F.33' + data$'F.34' + data$'F.44' + 
	                     (data$'G.13'/2) + (data$'G.14'/2) + (data$'G.23'/2) + (data$'G.24'/2) + data$'G.33' + data$'G.34' + data$'G.44' 
	data$diffMales  <-  (data$males - data$ZHat)

    # Color scheme
    COLS  <-  c(transparentColor('dodgerblue', opacity=0.8),
    	        transparentColor('dodgerblue4', opacity=0.8),
    	        transparentColor('tomato', opacity=0.8))

    # k index for easy plotting
    rs  <-  unique(data$r)
    ks  <-  unique(data$k)

    # Set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

##  Row 1: 
    # Panel 1: r = 0
    dat  <-  subset(data, data$r == rs[1])
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(min(data$diffMales),max(data$diffMales)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different 
        abline(h=0, lwd=2, col='black')
        lines(diffMales[k==ks[1]] ~ C[k==ks[1]], lwd=2, col=COLS[1], cex=1, data=dat)
        lines(diffMales[k==ks[2]] ~ C[k==ks[2]], lwd=2, col=COLS[2], cex=1, data=dat)
        lines(diffMales[k==ks[3]] ~ C[k==ks[3]], lwd=2, col=COLS[3], cex=1, data=dat)
#        lines(diffMales[k==ks[4]] ~ C[k==ks[4]], lwd=2, col=COLS[4], cex=1, data=dat)
        # axes
        axis(1, las=1,labels=NA)
        axis(2, las=1)
        proportionalLabel(0.5, 1.1, expression(paste(italic(r)," = 0.0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35, 0.5, expression(paste(Delta," Female frequency")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        #legend
        legend(
              x       =  usr[2]*0.975,
              y       =  usr[4],
              legend  =  c(
                          expression(paste(italic(k)~"="~italic(hat(k))%*%1.1)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%1.95)),
                          expression(paste(italic(k)~"="~italic(hat(k))%*%0.9))),
              lty     =  c(1,1,1),
              lwd     =  c(2,2,2),
              col     =  c(COLS[1],COLS[2],COLS[3]),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
    rm(dat)

    # Panel 2: r = 0.01
    dat  <-  subset(data, data$r == rs[2])
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(min(data$diffMales),max(data$diffMales)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different k values
        abline(h=0, lwd=2, col='black')
        lines(diffMales[k==ks[1]] ~ C[k==ks[1]], lwd=2, col=COLS[1], cex=1, data=dat)
        lines(diffMales[k==ks[2]] ~ C[k==ks[2]], lwd=2, col=COLS[2], cex=1, data=dat)
        lines(diffMales[k==ks[3]] ~ C[k==ks[3]], lwd=2, col=COLS[3], cex=1, data=dat)
        # axes
        axis(1, las=1,labels=NA)
        axis(2, las=1,labels=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(r)," = ", 0.02)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
    rm(dat)

##  Row 2: 
    # Panel 3: r = 0.05
    dat  <-  subset(data, data$r == rs[3])
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(min(data$diffMales),max(data$diffMales)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different 
        abline(h=0, lwd=2, col='black')
        lines(diffMales[k==ks[1]] ~ C[k==ks[1]], lwd=2, col=COLS[1], cex=1, data=dat)
        lines(diffMales[k==ks[2]] ~ C[k==ks[2]], lwd=2, col=COLS[2], cex=1, data=dat)
        lines(diffMales[k==ks[3]] ~ C[k==ks[3]], lwd=2, col=COLS[3], cex=1, data=dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.35, 0.5, expression(paste(Delta," Female frequency")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, 1.1, expression(paste(italic(r)," = ", 0.05)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.35, expression(paste(italic(C))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
    rm(dat)

    # Panel 4: r = 0.1
    dat  <-  subset(data, data$r == rs[4])
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.9), ylim = c(min(data$diffMales),max(data$diffMales)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Equilibrium frequencies for different 
        abline(h=0, lwd=2, col='black')
        lines(diffMales[k==ks[1]] ~ C[k==ks[1]], lwd=2, col=COLS[1], cex=1, data=dat)
        lines(diffMales[k==ks[2]] ~ C[k==ks[2]], lwd=2, col=COLS[2], cex=1, data=dat)
        lines(diffMales[k==ks[3]] ~ C[k==ks[3]], lwd=2, col=COLS[3], cex=1, data=dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(r)," = ", 0.1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.35, expression(paste(italic(C))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
    rm(dat)

}