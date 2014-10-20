'''
 Author: Kyle Cranmer   October 2014
 * Copyright (C) Kyle Cranmer
 * BSD License

This code quickly calculates the expected 95% CLs upper limit on the number of
signal events, s, given an expected background, bExp, and 
uncertainty on the background estimate, deltaB.
The background uncertainty is absolute (not relative) and is uncertainty
on the mean background (so you don't include Poisson fluctuatiosn in this number).
Example Usage: 
  * you expect 50 +/- 3 background events
  * ExpectedLimit(50,3)
  * returns s_95 = 13.4 events

The derivations of these formulae are based on a statistical model:
	Pois(n | s+b ) * Pois(m | tau * b)
The tau quantity is calculated from tau=bExp/deltaB/deltaB.
The maximum likelihood estimate and conditional maximum likelihood estimate
were solved analytically and coded here.
The log-likelihood ratio and the profile log likelihood ratio follow immediately.
The b-only and s+b p-values needed for CLs are calculated using the 
asymptotic distributions in Eur. Phys. J. C 71 (2011) 1554.  arXiv:1007.1727

'''

from __future__ import division
from scipy.stats import norm
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



def ExpectedLimit(bExp, deltaB) :

	tau = bExp/deltaB/deltaB

	s=0.1
	while CLs(bExp,bExp*tau,s,tau) < 0.95 :
		s+=0.01
	#print "Approximate expected 95% CLs upper-limit = ", s 

	# use http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html#scipy.optimize.brentq
	return s




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


