#####################################################
#  Reconsidering the evolution of dioecy from the
#  perspective of sexually antagonistic selection
#
#  Necessary functions for analyses of the model
#  of the genotypic frequency recursions
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
## Change in genotypic frequencies due to fitness
##  through female function
Ff.11  <-  function(Fii, Wf.mat, ...) {
	Fii[1]*Wf.mat[1,1]
}
Ff.12  <-  function(Fii, Wf.mat, ...) {
	Fii[2]*Wf.mat[1,2]
}
Ff.13  <-  function(Fii, Wf.mat, ...) {
	Fii[3]*Wf.mat[1,3]
}
Ff.14  <-  function(Fii, Wf.mat, ...) {
	Fii[4]*Wf.mat[1,4]
}
Ff.22  <-  function(Fii, Wf.mat, ...) {
	Fii[5]*Wf.mat[2,2]
}
Ff.23  <-  function(Fii, Wf.mat, ...) {
	Fii[6]*Wf.mat[2,3]
}
Ff.24  <-  function(Fii, Wf.mat, ...) {
	Fii[7]*Wf.mat[2,4]
}
Ff.33  <-  function(Fii, Wf.mat, ...) {
	Fii[8]*Wf.mat[3,3]
}
Ff.34  <-  function(Fii, Wf.mat, ...) {
	Fii[9]*Wf.mat[3,4]
}
Ff.44  <-  function(Fii, Wf.mat, ...) {
	Fii[10]*Wf.mat[4,4]
}


#######################################################
## Change in genotypic frequencies due to fitness
##  through male function
Fm.11  <-  function(Fii, Wm.mat, ...) {
	Fii[1]*Wm.mat[1,1]
}
Fm.12  <-  function(Fii, Wm.mat, ...) {
	Fii[2]*Wm.mat[1,2]
}
Fm.13  <-  function(Fii, Wm.mat, ...) {
	Fii[3]*Wm.mat[1,3]
}
Fm.14  <-  function(Fii, Wm.mat, ...) {
	Fii[4]*Wm.mat[1,4]
}
Fm.22  <-  function(Fii, Wm.mat, ...) {
	Fii[5]*Wm.mat[2,2]
}
Fm.23  <-  function(Fii, Wm.mat, ...) {
	Fii[6]*Wm.mat[2,3]
}
Fm.24  <-  function(Fii, Wm.mat, ...) {
	Fii[7]*Wm.mat[2,4]
}
Fm.33  <-  function(Fii, Wm.mat, ...) {
	Fii[8]*Wm.mat[3,3]
}
Fm.34  <-  function(Fii, Wm.mat, ...) {
	Fii[9]*Wm.mat[3,4]
}
Fm.44  <-  function(Fii, Wm.mat, ...) {
	Fii[10]*Wm.mat[4,4]
}

#######################################################
## Change in genotypic frequencies due to fitness
##  through selfing
Fs.11  <-  function(Fii, Ws.mat, ...) {
	Fii[1]*Ws.mat[1,1]
}
Fs.12  <-  function(Fii, Ws.mat, ...) {
	Fii[2]*Ws.mat[1,2]
}
Fs.13  <-  function(Fii, Ws.mat, ...) {
	Fii[3]*Ws.mat[1,3]
}
Fs.14  <-  function(Fii, Ws.mat, ...) {
	Fii[4]*Ws.mat[1,4]
}
Fs.22  <-  function(Fii, Ws.mat, ...) {
	Fii[5]*Ws.mat[2,2]
}
Fs.23  <-  function(Fii, Ws.mat, ...) {
	Fii[6]*Ws.mat[2,3]
}
Fs.24  <-  function(Fii, Ws.mat, ...) {
	Fii[7]*Ws.mat[2,4]
}
Fs.33  <-  function(Fii, Ws.mat, ...) {
	Fii[8]*Ws.mat[3,3]
}
Fs.34  <-  function(Fii, Ws.mat, ...) {
	Fii[9]*Ws.mat[3,4]
}
Fs.44  <-  function(Fii, Ws.mat, ...) {
	Fii[10]*Ws.mat[4,4]
}

#####################################################
## Average fitness through each sex role and selfing
Wf.av  <-  function(Fii, Wf.mat, ...){
   Ff.11(Fii,Wf.mat) + Ff.12(Fii,Wf.mat)+ Ff.13(Fii,Wf.mat) + Ff.14(Fii,Wf.mat) + Ff.22(Fii,Wf.mat) + Ff.23(Fii,Wf.mat) + Ff.24(Fii,Wf.mat) + Ff.33(Fii,Wf.mat) + Ff.34(Fii,Wf.mat) + Ff.44(Fii,Wf.mat)
}
Wm.av  <-  function(Fii, Wm.mat, ...){
   Fm.11(Fii,Wm.mat) + Fm.12(Fii,Wm.mat)+ Fm.13(Fii,Wm.mat) + Fm.14(Fii,Wm.mat) + Fm.22(Fii,Wm.mat) + Fm.23(Fii,Wm.mat) + Fm.24(Fii,Wm.mat) + Fm.33(Fii,Wm.mat) + Fm.34(Fii,Wm.mat) + Fm.44(Fii,Wm.mat)
}
Ws.av  <-  function(Fii, Ws.mat, ...){
   Fs.11(Fii,Ws.mat) + Fs.12(Fii,Ws.mat)+ Fs.13(Fii,Ws.mat) + Fs.14(Fii,Ws.mat) + Fs.22(Fii,Ws.mat) + Fs.23(Fii,Ws.mat) + Fs.24(Fii,Ws.mat) + Fs.33(Fii,Ws.mat) + Fs.34(Fii,Ws.mat) + Fs.44(Fii,Ws.mat)
}


#########################################
## Haplotype frequency change in gametes

# Ovules
x1  <-  function(Fii, Wf.mat, par.list, ...) {
	((2*Ff.11(Fii,Wf.mat) + Ff.12(Fii,Wf.mat) + Ff.13(Fii,Wf.mat) + Ff.14(Fii,Wf.mat)) / (2*Wf.av(Fii, Wf.mat))) - 
		par.list$r*((Ff.14(Fii,Wf.mat) - Ff.23(Fii,Wf.mat))/(2*Wf.av(Fii, Wf.mat)))
}
x2  <-  function(Fii, Wf.mat, par.list, ...){
	((2*Ff.22(Fii,Wf.mat) + Ff.12(Fii,Wf.mat) + Ff.23(Fii,Wf.mat) + Ff.24(Fii,Wf.mat))/(2*Wf.av(Fii, Wf.mat))) + 
		par.list$r*((Ff.14(Fii,Wf.mat) - Ff.23(Fii,Wf.mat))/(2*Wf.av(Fii, Wf.mat)))
}
x3  <-  function(Fii, Wf.mat, par.list, ...) {
	 ((2*Ff.33(Fii,Wf.mat) + Ff.34(Fii,Wf.mat) + Ff.13(Fii,Wf.mat) + Ff.23(Fii,Wf.mat))/(2*Wf.av(Fii, Wf.mat))) + 
 		par.list$r*((Ff.14(Fii,Wf.mat) - Ff.23(Fii,Wf.mat))/(2*Wf.av(Fii, Wf.mat)))
}
x4  <-  function(Fii, Wf.mat, par.list, ...) {
	((2*Ff.44(Fii,Wf.mat) + Ff.34(Fii,Wf.mat) + Ff.14(Fii,Wf.mat) + Ff.24(Fii,Wf.mat))/(2*Wf.av(Fii, Wf.mat))) - 
		par.list$r*((Ff.14(Fii,Wf.mat) - Ff.23(Fii,Wf.mat))/(2*Wf.av(Fii, Wf.mat)))
}


# Pollen/Sperm
y1  <-  function(Fii, Wm.mat, par.list, ...) {
	((2*Fm.11(Fii,Wm.mat) + Fm.12(Fii,Wm.mat) + Fm.13(Fii,Wm.mat) + Fm.14(Fii,Wm.mat)) / (2*Wm.av(Fii, Wm.mat))) - 
		par.list$r*((Fm.14(Fii,Wm.mat) - Fm.23(Fii,Wm.mat))/(2*Wm.av(Fii, Wm.mat)))
}
y2  <-  function(Fii, Wm.mat, par.list, ...){
	((2*Fm.22(Fii,Wm.mat) + Fm.12(Fii,Wm.mat) + Fm.23(Fii,Wm.mat) + Fm.24(Fii,Wm.mat))/(2*Wm.av(Fii, Wm.mat))) + 
		par.list$r*((Fm.14(Fii,Wm.mat) - Fm.23(Fii,Wm.mat))/(2*Wm.av(Fii, Wm.mat)))
}
y3  <-  function(Fii, Wm.mat, par.list, ...) {
	 ((2*Fm.33(Fii,Wm.mat) + Fm.34(Fii,Wm.mat) + Fm.13(Fii,Wm.mat) + Fm.23(Fii,Wm.mat))/(2*Wm.av(Fii, Wm.mat))) + 
 		par.list$r*((Fm.14(Fii,Wm.mat) - Fm.23(Fii,Wm.mat))/(2*Wm.av(Fii, Wm.mat)))
}
y4  <-  function(Fii, Wm.mat, par.list, ...) {
	((2*Fm.44(Fii,Wm.mat) + Fm.34(Fii,Wm.mat) + Fm.14(Fii,Wm.mat) + Fm.24(Fii,Wm.mat))/(2*Wm.av(Fii, Wm.mat))) - 
		par.list$r*((Fm.14(Fii,Wm.mat) - Fm.23(Fii,Wm.mat))/(2*Wm.av(Fii, Wm.mat)))
}




#########################################
## Genotypic frequency recursions
F11.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list,...) {
	x1  <-  x1(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y1  <-  y1(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C) * x1*y1 + 
	     par.list$C  * ((Fs.11(Fii,Ws.mat) + Fs.12(Fii,Ws.mat)/4 + Fs.13(Fii,Ws.mat)/4 + Fs.14(Fii,Ws.mat)*((1 - par.list$r)^2)/4 + Fs.23(Fii,Ws.mat)*(par.list$r^2)/4)/(Ws.av(Fii=Fii, Ws.mat=Ws.mat)))
}
F12.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list,...) {
	x1  <-  x1(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y2  <-  y2(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	x2  <-  x2(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y1  <-  y1(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)* (x1*y2 + x2*y1) + 
	     par.list$C * ((Fs.12(Fii,Ws.mat)/2 + Fs.14(Fii,Ws.mat)*par.list$r*(1 - par.list$r)/2 + Fs.23(Fii,Ws.mat)*par.list$r*(1 - par.list$r)/2)/(Ws.av(Fii=Fii, Ws.mat=Ws.mat)))
}
F13.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list,...) {
	x1  <-  x1(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y3  <-  y3(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	x3  <-  x3(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y1  <-  y1(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C) * (x1*y3 + x3*y1) + 
	     par.list$C  * ((Fs.13(Fii,Ws.mat)/2 + Fs.14(Fii,Ws.mat)*par.list$r*(1 - par.list$r)/2 + Fs.23(Fii,Ws.mat)*par.list$r*(1 - par.list$r)/2)/(Ws.av(Fii=Fii, Ws.mat=Ws.mat)))
}
F14.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list,...) {
	x1  <-  x1(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y4  <-  y4(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	x4  <-  x4(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y1  <-  y1(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C) * (x1*y4 + x4*y1) + 
	     par.list$C  * ((Fs.14(Fii,Ws.mat)*((1 - par.list$r)^2)/2 + Fs.23(Fii,Ws.mat)*(par.list$r^2)/2)/(Ws.av(Fii=Fii, Ws.mat=Ws.mat)))
}
F22.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list,...) {
	x2  <-  x2(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y2  <-  y2(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C) * x2*y2 + 
	     par.list$C  * ((Fs.22(Fii,Ws.mat) + Fs.12(Fii,Ws.mat)/4 + Fs.14(Fii,Ws.mat)*(par.list$r^2)/4 + Fs.23(Fii,Ws.mat)*((1 - par.list$r)^2)/4 + Fs.24(Fii,Ws.mat)/4)/(Ws.av(Fii=Fii, Ws.mat=Ws.mat)))
}
F23.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list,...) {
	x2  <-  x2(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y3  <-  y3(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	x3  <-  x3(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y2  <-  y2(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C) * (x2*y3 + x3*y2) + 
	     par.list$C  * ((Fs.14(Fii,Ws.mat)*(par.list$r^2)/2 + Fs.23(Fii,Ws.mat)*((1 - par.list$r)^2)/2)/(Ws.av(Fii=Fii, Ws.mat=Ws.mat)))
}
F24.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list,...) {
	x2  <-  x2(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y4  <-  y4(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	x4  <-  x4(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y2  <-  y2(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C) * (x2*y4 + x4*y2) + 
	     par.list$C  * ((Fs.24(Fii,Ws.mat)/2 + Fs.14(Fii,Ws.mat)*par.list$r*(1 - par.list$r)/2 + Fs.23(Fii,Ws.mat)*par.list$r*(1 - par.list$r)/2)/(Ws.av(Fii=Fii, Ws.mat=Ws.mat)))
}
F33.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list,...) {
	x3  <-  x3(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y3  <-  y3(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*x3*y3 + par.list$C*((Fii[8]*Wf.mat[3,3] + Fii[3]*Wf.mat[1,3]/4 + Fii[4]*Wf.mat[1,4]*(par.list$r^2)/4 + Fii[6]*Wf.mat[2,3]*((1 - par.list$r)^2)/4 + Fii[9]*Wf.mat[3,4]/4)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}
F34.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list,...) {
	x3  <-  x3(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y4  <-  y4(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	x4  <-  x4(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y3  <-  y3(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C) * (x3*y4 + x4*y3) + 
	     par.list$C  * ((Fs.14(Fii,Ws.mat)*par.list$r*(1 - par.list$r)/2 + Fs.23(Fii,Ws.mat)*par.list$r*(1 - par.list$r)/2 + Fs.34(Fii,Ws.mat)/2)/(Ws.av(Fii=Fii, Ws.mat=Ws.mat)))
}
F44.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list,...) {
	x4  <-  x4(Fii = Fii, Wf.mat = Wf.mat, par.list = par.list,...)
	y4  <-  y4(Fii = Fii, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C) * x4*y4 + 
	     par.list$C  * ((Fs.44(Fii,Ws.mat) + Fs.14(Fii,Ws.mat)*((1 - par.list$r)^2)/4 + Fs.23(Fii,Ws.mat)*(par.list$r^2)/4 + Fs.24(Fii,Ws.mat)/4 + Fs.34(Fii,Ws.mat)/4)/(Ws.av(Fii=Fii, Ws.mat=Ws.mat)))
}


################################################
## Single-locus equilibrium frequency functions

qHatAdd  <-  function(C, sf, sm) {
	(sf*(1+C)-sm*(1-C)*(1-sf))/(2*sf*sm)
}

#qHatDomRev  <-  function(C, sf, sm, hf, hm) {
#}

QEFAA  <-  function(q,C) {
	(1-q)^2 + (C*q*(1-q))/(2-C)
}
QEFAa  <-  function(q,C) {
	2*q*(1-q) - (C*q*(1-q))/(2-C)
}
QEFaa  <-  function(q,C) {
	q^2 + (C*q*(1-q))/(2-C)
}

# Kidwell et al. (1977) male and female frequencies
pmHat  <-  function(sf,sm) {
	((sm-1)/sm) + sqrt((sm*sf-sm-sf+2)/(2*sm*sf))
}

pfHat  <-  function(sf,sm) {
	-(-2+sf * sqrt((4+2*sf*(sm-1)-2*sm)/(sf*sm)))/(2*sf)
}
KidwellFAA  <-  function(pf,pm) {
	pf*pm
}
KidwellFAa  <-  function(pf,pm) {
	pf*(1-pm) + (1-pf)*pm
}
KidwellFaa  <-  function(pf,pm) {
	(1-pf)*(1-pm)
}













########################
##  Simulation function
########################

#' Forward deterministic simulation of genotypic recursions for
#' 2-locus SA w/ linked sterility locus
#'
#' @title Forward deterministic simulation of genotypic recursions
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'				   gen  =  5000,
#'				   C    =  0,
#'				   sm   =  0.7,
#'				   sf   =  0.7,
#'				   hm   =  1/2,
#'				   hf   =  1/2,
#'                 hMf  =  0,
#'                 hMm  =  0,
#'                 k    =  1,
#'				   r    =  0.5
#'				   )
#' @param Fii.init A vector of initial genotypic frequencies (must have length = 10).
#' c(0.99,0,0,0,0,0,0,0,0,0.01) for invasion of aabb into population 'fixed' for AABB.
#' c(0.01,0,0,0,0,0,0,0,0,0.99) for invasion of AABB into population 'fixed' for aabb.
#' @return Returns a list with timeseries for each genotype, equilibrium frequencies, and a numeric (0,1) for whether the 
#' equilibrium was polymorphic (with tolerance 1E-6).
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSim(par.list, Fii.init, threshold = 1e-6) 
recursionFwdSim  <-  function(par.list, Fii.init, threshold = 1e-6) {

	##  Warnings
	if(any(par.list[2:8] < 0) | any(par.list[2:8] > 1) | par.list[10] > 0.5)
		stop('The chosen parameter values fall outside of the reasonable bounds')

	if(par.list$hf  !=  par.list$hm)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(par.list$hf != 0.5 & par.list$hf != 0.25)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')


	##  Fitness Matrices
	Wf.mat   <-  matrix(
						c(1,                                                    (1 + par.list$k*par.list$hMf),                         (1-par.list$hf*par.list$sf),                          (1-par.list$hf*par.list$sf + par.list$k*par.list$hMf),
						 (1 + par.list$k*par.list$hMf),                         (1 + par.list$k),                                      (1-par.list$hf*par.list$sf + par.list$k*par.list$hMf), (1-par.list$hf*par.list$sf + par.list$k),
						 (1-par.list$hf*par.list$sf),                           (1-par.list$hf*par.list$sf + par.list$k*par.list$hMf), (1-par.list$sf),                                       (1-par.list$sf + par.list$k*par.list$hMf),
						 (1-par.list$hf*par.list$sf + par.list$k*par.list$hMf), (1-par.list$hf*par.list$sf + par.list$k),              (1-par.list$sf + par.list$k*par.list$hMf),            (1-par.list$sf + par.list$k)), 
						nrow=4, byrow=TRUE
						)

	Wm.mat   <-  matrix(
						c((1-par.list$sm),                              (1-par.list$sm)*(1-par.list$hMm),             (1-par.list$hm*par.list$sm),                  (1-par.list$hm*par.list$sm)*(1-par.list$hMm),
						  (1-par.list$sm)*(1-par.list$hMm),             0,                                            (1-par.list$hm*par.list$sm)*(1-par.list$hMm), 0,
						  (1-par.list$hm*par.list$sm),                  (1-par.list$hm*par.list$sm)*(1-par.list$hMm), 1,                                            (1-par.list$hMm),
						  (1-par.list$hm*par.list$sm)*(1-par.list$hMm), 0,                                            (1-par.list$hMm),                             0), 
						nrow=4, byrow=TRUE
						)	

	Ws.mat   <-  matrix(
						c(1,                                                                     (1 + par.list$k*par.list$hMf)*(1-par.list$hMm),                         (1-par.list$hf*par.list$sf) ,                                           (1-par.list$hf*par.list$sf + par.list$k*par.list$hMf)*(1-par.list$hMm),
						 (1 + par.list$k*par.list$hMf)*(1-par.list$hMm),                         0,                                                                      (1-par.list$hf*par.list$sf + par.list$k*par.list$hMf)*(1-par.list$hMm), 0,
						 (1-par.list$hf*par.list$sf),                                            (1-par.list$hf*par.list$sf + par.list$k*par.list$hMf)*(1-par.list$hMm), (1-par.list$sf),                                                        (1-par.list$sf + par.list$k*par.list$hMf)*(1-par.list$hMm),
						 (1-par.list$hf*par.list$sf + par.list$k*par.list$hMf)*(1-par.list$hMm), 0,                                                                      (1-par.list$sf + par.list$k*par.list$hMf)*(1-par.list$hMm),            0), 
						nrow=4, byrow=TRUE
						)


	##  Initilize data storage structures
	Fii.gen  <-  matrix(0, ncol=10, nrow=par.list$gen)

	##  Initial frequencies
#	if(par.list$sf > par.list$sm)
#		Fii.init    <-  c(0.01,0,0,0,0,0,0,0,0,0.99)
#	if(par.list$sf < par.list$sm)
#		Fii.init    <-  c(0.99,0,0,0,0,0,0,0,0,0.01)
    

	##  Generation Loop
		# initialize
		Fii.gen[1,1]   <-  round(F11.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[1,2]   <-  round(F12.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[1,3]   <-  round(F13.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[1,4]   <-  round(F14.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[1,5]   <-  round(F22.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[1,6]   <-  round(F23.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[1,7]   <-  round(F24.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[1,8]   <-  round(F33.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[1,9]   <-  round(F34.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[1,10]  <-  round(F44.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)


	# Start simulation
	i      <-  2
	diffs  <-  rep(1,10)

	while (i < par.list$gen & any(diffs[diffs != 0] > threshold)) {
		Fii.gen[i,1]   <-  round(F11.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[i,2]   <-  round(F12.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[i,3]   <-  round(F13.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[i,4]   <-  round(F14.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[i,5]   <-  round(F22.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[i,6]   <-  round(F23.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[i,7]   <-  round(F24.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[i,8]   <-  round(F33.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[i,9]   <-  round(F34.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		Fii.gen[i,10]  <-  round(F44.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, Ws.mat = Ws.mat, par.list = par.list), digits=6)
		
		diffs  <-  Fii.gen[i,] - Fii.gen[i-1,]
		i      <-  i+1
	}

	##  Is equilibrium polymorphic?
#	if (any(Fii.gen[i-1,c(1,10)] > 0.999)) # & all(Fii.gen[i-1,2:9] < 1e-4))
#		 Poly  <-  0
#	else Poly  <-  1

	##  Calculate Eigenvalues from analytic solutions 
	##  using quasi-equibirium genotypic frequencies
#	if (par.list$hf == 0.5) {
#		l.AB1  <- lambda.AB1.add(par.list)
#		l.AB2  <- lambda.AB2.add(par.list)
#		l.ab1  <- lambda.ab1.add(par.list)
#		l.ab2  <- lambda.ab2.add(par.list)
#	}
#	if (par.list$hf == 0.25) {
#		l.AB1  <- lambda.AB1.domRev(par.list)
#		l.AB2  <- lambda.AB2.domRev(par.list)
#		l.ab1  <- lambda.ab1.domRev(par.list)
#		l.ab2  <- lambda.ab2.domRev(par.list)
#	}

#	if (any(c(l.AB1, l.AB2) > 1) & any(c(l.ab1, l.ab2) > 1 ))
#		 eigPoly  <-  1
#	else eigPoly  <-  0

	##  Does simulation result agree with Eigenvalues?

#	if (Poly == eigPoly)
#		 agree  <-  1
#	else agree  <-  0


	##  Output list
	res  <-  list(
				  "par.list" =  par.list,
				  "Fii.gen"  =  Fii.gen[1:i-1,],
				  "EQ.freq"  =  Fii.gen[i-1,]
#				  "l.AB1"    =  l.AB1,
#				  "l.AB2"    =  l.AB2,
#				  "l.ab1"    =  l.ab1,
#				  "l.ab2"    =  l.ab2,
#				  "Poly"     =  Poly,
#				  "eigPoly"  =  eigPoly,
#				  "agree"    =  agree
 				 )
	return(res)
}

