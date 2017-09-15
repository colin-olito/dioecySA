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



###############################################################
## General Functions 
##    -- Single-locus SA equilibrium allele frequencies
##    -- Single-locus SA equilibrium genotypic frequencies
##    -- Single-locus SA Invasion Conditions
##    -- Single-locus Sterility allele equilibrium frequencies
###############################################################

# Single-locus SA equilibrium allele frequencies
qHatAdd  <-  function(C, delta, sf, sm) {
#	qHat  <-  ((sf - sm + sf*sm) + C*(sf + sm - sf*sm - 2*sf*delta)) / (2*(sf*sm - C*sf*sm*delta))
	qHat  <-  (sf - sm + sf*sm + C*(sf + sm - sf*sm - 2*sf*delta)) /(2*(sf*sm - C*sf*sm*delta))
# 	if(qHat < 0 | qHat > 1 | is.nan(qHat) | is.na(qHat))
#		qHat  <-  0
	qHat
}
qHatDomRev  <-  function(C, delta, sf, sm, h) {
	qHat  <-  ((-1 + C)*sm*(-2*h + C*(-1 + 2*h + delta)) + sf*(2 - 2*h + C*(-1 + 2*h - delta))*(-1 + C*(-1 + 2*delta))) / (2*(-1 + C)*(-1 + 2*h)*((-1 + C)*sm + sf*(-1 + C*(-1 + 2*delta))))
	if(qHat < 0 | qHat > 1 | is.nan(qHat) | is.na(qHat))
		qHat  <-  0
	qHat
}

# Single-locus SA equilibrium genotypic frequencies
QE.FAA  <-  function(q, C) {
	(1-q)^2 + (C*q*(1-q))/(2-C)
}
QE.FAa  <-  function(q,C) {
	2*q*(1-q) - (C*q*(1-q))/(2-C)
}
QE.Faa  <-  function(q,C) {
	q^2 + (C*q*(1-q))/(2-C)
}

# Kidwell et al. (1977) male and female frequencies
pmHat  <-  function(sf,sm) {
	pmHat  <-  ((sm-1)/sm) + sqrt((sm*sf-sm-sf+2)/(2*sm*sf))
	if(pmHat < 0 | pmHat > 1 | is.nan(pmHat) | is.na(pmHat))
		pmHat  <-  0
	pmHat
}
pfHat  <-  function(sf,sm) {
	pfHat  <-  -(-2+sf * sqrt((4+2*sf*(sm-1)-2*sm)/(sf*sm)))/(2*sf)
	if(pfHat < 0 | pfHat > 1 | is.nan(pfHat) | is.na(pfHat))
		pfHat  <-  0
	pfHat
}
Kidwell.FAA  <-  function(pf,pm) {
	pf*pm
}
Kidwell.FAa  <-  function(pf,pm) {
	pf*(1-pm) + (1-pf)*pm
}
Kidwell.Faa  <-  function(pf,pm) {
	(1-pf)*(1-pm)
}

## Single-locus SA invasion conditions
Inv.a.add  <-  function(sm, C, delta) {
	((1 - C)*sm) / ((sm - 1)*(-1 + C*(-1 + 2*delta)))
}
Inv.A.add  <-  function(sm, C, delta) {
	(sm - C*sm) / (1 + C + sm - C*sm - 2*C*delta)
}
Inv.a.domRev  <-  function(sm, h, C, delta) {
	((-1 + C)*(2 - C + 2*(-1 + C)*h)*sm) / ((-C + 2*(-1 + C)*h)*(-1 + sm)*(-1 + C*(-1 + 2*delta)))
}
Inv.A.domRev  <-  function(sm, h, C, delta) {
	((-1 + C)*(-C + 2*(-1 + C)*h)*sm) / (2 + 2*h*(-1 + sm) + (C^2)*(-1 + 2*h)*(1 + sm - 2*delta) + C*(1 + sm - 4*h*sm + 4*(-1 + h)*delta))
}
Inv.a.compDomRev  <-  function(sm, C, delta) {
	((C - 2)*(C - 1)*sm)/(C*(sm - 1)*(C*(2*delta - 1) - 1))
}
Inv.A.compDomRev  <-  function(sm, C, delta) {
	((C - 1)*C*sm)/(-2 + (C - 1)*C*(sm + 1) - 2*(C - 2)*C*delta)
}

#############################################################
# Relevant equations from Charlesworth & Charlesworth (1978)
# Single-locus sterility allele equilibrium frequencies, Eqs.5 & 9 in Charlesworth & Charlesworth (1978).
# Designed for use with recursion simulations, so takes par.list as argument rather than separate args for
# each parameter k, C, delta.
Zhat.gyn  <-  function(k, C, delta) {
	Zhat  <-  (k + 2*C*delta - 1) / (2*(k + C*delta))
	if(Zhat < 0) {
		Zhat  <-  0
	}
	Zhat
}
Zhat.gyn.list  <-  function(par.list) {
	Zhat  <-  (par.list$k + 2*par.list$C*par.list$delta - 1) / (2*(par.list$k + par.list$C*par.list$delta))
	if(Zhat < 0) {
		Zhat  <-  0
	}
	Zhat
}

Zhat.and  <-  function(k, C, delta) {
	Zhat  <-  ((1 + k)*(1 - C) - 2*(1 - C*delta)) / (2*k*(1 - C*delta))
	if(Zhat < 0){
		Zhat  <-  0
	}
	Zhat
}
Zhat.and.list  <-  function(par.list) {
	Zhat  <-  ((1 + par.list$k)*(1 - par.list$C) - 2*(1 - par.list$C*par.list$delta)) / (2*par.list$k*(1 - par.list$C*par.list$delta))
	if(Zhat < 0){
		Zhat  <-  0
	}
	Zhat
}

# Single-locus invasion criteria for sterility alleless (Gynodioecy & Androdioecy).
# Correspond to Eqs.4 & 8 in Charlesworth * Charlesworth (1978)
invGyn  <-  function(C, delta) {
	1 - 2*C*delta
}

invAnd  <-  function(C, delta) {
	(1 + C - 2*C*delta) / (1 - C)
}

#############################################################
# Functions to account for C*delta correlation for recursion 
# simulation figures.

L.phenom  <- function(C, a) {
	((a*(1 - C))/(C + a*(1 - C)))
}

deltaC  <-  function(dStar = 1, C = 0, a = 1, b = 0.5) {
	dStar*(1 - b*(1 - L.phenom(C,a)))
}


#' Weak selection approximation for equal frequencies of A and a alleles at SA locus under additive fitness effects
#' 
#' @title Weak selection approximation for equal frequencies of A and a alleles at SA locus 
#' @param sm Selection coefficient through male sex function
#' @param C Population selfing rate
#' @param delta Population inbreeding depression
#' @return Returns a single value for the corresponding value of sf
#' @seealso `qHatAdd`
#' @export
#' @author Colin Olito.
#'
equalPQ.add  <-  function(sm, C, delta) {
	(sm - C*sm)/(1 + C - C*sm - 2*C*delta + C*sm*delta)
}

#' Weak selection approximation for equal frequencies of A and a alleles at SA locus under dominance reversal fitness effects
#' 
#' @title Weak selection approximation for equal frequencies of A and a alleles at SA locus 
#' @param sm Selection coefficient through male sex function
#' @param C Population selfing rate
#' @param delta Population inbreeding depression
#' @return Returns a single value for the corresponding value of sf
#' @seealso `qHatDomRev`
#' @export
#' @author Colin Olito.
#'
equalPQ.domRev  <-  function(sm, C, delta) {
	((C - 1)*sm)/(-1 - C + 2*C*delta)
}
