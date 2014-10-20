'''
 Author: Kyle Cranmer   October 2014
 * Copyright (C) Kyle Cranmer
_ License LGPL v2.1 (for ROOT compatibility, happy to make it BSD for other purposes) _

This code quickly calculates the expected 95% CLs upper limit on the number of
signal events, s, given an expected background, bExp, and 
uncertainty on the background estimate, deltaB.
The background uncertainty is absolute (not relative) and is uncertainty
on the mean background (so you don't include sqrt(bExp) Poisson fluctuatiosn in this number).

Example Usage: 
  * you expect 50 +/- 3 background events
  * ExpectedLimit(50,3)
  * returns s_95 = 13.4 events

Similar code for expected discovery significance is also included
Example Usage:
	* you expect 50 signal events, 100 +/- 7 background events
	* ExpectedSignificance(50,100,7)
	* returns 3.72 sigma

The derivations of these formulae are based on a statistical model:
	Pois(n | s+b ) * Pois(m | tau * b)
The tau quantity is calculated from tau=bExp/deltaB/deltaB.
See arXiv:physics/0702156 for motivation of this model.
The maximum likelihood estimate and conditional maximum likelihood estimate
were solved analytically and coded here.
The log-likelihood ratio and the profile log likelihood ratio follow immediately.
The b-only and s+b p-values needed for CLs are calculated using the 
asymptotic distributions in 
	Cowan, Cranmer, Gross, Vitells,
	Eur. Phys. J. C 71 (2011) 1554.  
	arXiv:1007.1727.
The ExpectedSignificance formulae was derived by Cowan and the numerical solution 
for the Expected upper limit was written by Cranmer as correlaries to that paper.

Note, the ExpectedSignificance euqation is the same as Eq.(17) of
 Tipei Li and Yuqian Ma, Astrophysical Journal 272 (1983) 317–324.
and Eq.(25) of 
Robert D. Cousins, James T. Linnemann and Jordan Tucker, NIM A 595 (2008) 480– 501; arXiv:physics/0702156.
after making the replacements n=bExp and m=bExp*tau. 
'''

from __future__ import division
from scipy.stats import norm
from scipy.optimize import brentq
from numpy import sqrt, log

def bhathat(n, m, s, tau):
	return (m+n-s-s*tau + sqrt(4*m*s*(1+tau) + pow(s+s*tau-m-n,2) ) ) / (2 * (1+tau )) 


def shat(n, m, tau):
	return (n*tau -m )/tau


def bhat(n, m, tau):
	return m/tau


def sigma(n, m, s, tau):
	return 8*m*n*pow(1+tau,2)* \
		(m*m+2*m*n+n*n+m*s-n*s+m*s*tau-n*s*tau+(n+m)*sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2))) \
		/sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)) \
		/pow(m+n-s-s*tau+sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)),2) \
		/pow(m+n+s+s*tau+sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)),2) 


def logL(n, m, s, b, tau):
	return (s+b) - n * log(s+b) + (tau*b) - m * log(tau*b)


def logLambda(n, m, s, tau):
	return logL(n, m, s, bhathat(n,m,s,tau), tau)  \
		- logL(n,m,shat(n,m,tau), bhat(n,m,tau), tau)


def F(qmu, mu, muprime, sigma):
	return norm.cdf(sqrt(qmu)-(mu-muprime)/sigma)


def CLs(n, m, s, tau):
	qmu = 2*logLambda(n,m,s,tau)
	sig = sigma(n,m,s,tau)
	if shat(n,m,tau) > s :
		return 0.5

	return F(qmu,s,s,sig)/(1-F(qmu,s,0,sig))


#
'''
class CLsHelper() : 
	n=0
	m=0
	tau=0
	evalCLs (s ) const : 
   		return CLs(n,m,s,tau)
	
'''

def CLsArgumentWrapper(s,n,m,tau):
	return CLs(n,m,s,tau)-0.95


def ExpectedLimit(bExp, deltaB) :

	tau = bExp/deltaB/deltaB

	xhi = 3*sqrt(bExp)+3*deltaB #a reasonable guess of what is larger than the upper limit
	try:
		s = brentq(CLsArgumentWrapper, 0.001, xhi, args=(bExp,bExp*tau,tau))
		print "Approximate expected 95% CLs upper-limit = ", s 
	except:
		print "brentq failed (boundaries?) using a simple scan"
		s=0.01
		while CLs(bExp,bExp*tau,s,tau) < 0.95 :
			s+=0.1
		print "Approximate expected 95% CLs upper-limit = ", s 
		return s
 

	# use http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html#scipy.optimize.brentq




def ExpectedSignificance(s, b, deltaB):
	# same idea, but this is expected significance (in sigma)
	return sqrt(2*( (s+b)*log( ((s+b)*(b+deltaB*deltaB)) / (b*b+(s+b)*deltaB*deltaB) ) \
				-( b*b/deltaB/deltaB)*log(1.+ (s*deltaB*deltaB)/(b*(b+deltaB*deltaB)) ) 
				  )
			)


if __name__ == '__main__':

	if	shat(100,50,1)==50 :
		print  "ok" 
	else:
		print  "oops" 

	if	bhat(100,50,1)==50 :
		print  "ok" 
	else:
		print  "oops" 

	if 	bhathat(100,50,50,1)==50 :
		print  "ok" 
	else:
		print  "oops" 

	if	149 < 1/sigma(100,50,50,1) and 1/sigma(100,50,50,1) <151 :
		print  "ok" 
	else:
		print  "oops" 

	if	logLambda(100,50,50,1)==0 :
		print  "ok" 
	else:
		print  "oops" 

	if	ExpectedSignificance(50,100,10)>3. and ExpectedSignificance(50,100,10)<3.2 :
		print  "ok" 
	else:
		print  "oops" 

	if	ExpectedLimit(50,7)>16.4 and ExpectedLimit(50,7)<16.5 :
		print  "ok" 
	else:
		print  "oops" 


