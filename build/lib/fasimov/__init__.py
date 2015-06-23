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
  * you expect 50 +/- 7 background events
  * ExpectedLimit(50,7)
  * returns s_95 = 19.7 events

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
Tipei Li and Yuqian Ma, Astrophysical Journal 272 (1983) 317-324
and Eq.(25) of
Robert D. Cousins, James T. Linnemann and Jordan Tucker, NIM A 595 (2008) 480-501; arXiv:physics/0702156.
after making the replacements n=bExp and m=bExp*tau. 
'''

from __future__ import division
from scipy.stats import norm
from scipy.optimize import brentq
from numpy import sqrt, log, abs

def bhathat(n, m, s, tau):
	return (m+n-s-s*tau + sqrt(4*m*s*(1+tau) + pow(s+s*tau-m-n,2) ) ) / (2 * (1+tau )) 


def shat(n, m, tau):
	return (n*tau -m )/tau


def bhat(n, m, tau):
	return m/tau

def logL(n, m, s, b, tau):
	return (s+b) - n * log(s+b) + (tau*b) - m * log(tau*b)


def logLambda(n, m, s, tau):
	return logL(n, m, s, bhathat(n,m,s,tau), tau)  \
		- logL(n,m,shat(n,m,tau), bhat(n,m,tau), tau)

def fisher(n,m,s,tau):
	#from fisher information matrix approach
	return 8*m*n*pow(1+tau,2)* \
		(m*m+2*m*n+n*n+m*s-n*s+m*s*tau-n*s*tau+(n+m)*sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2))) \
		/sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)) \
		/pow(m+n-s-s*tau+sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)),2) \
		/pow(m+n+s+s*tau+sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)),2) 

def sigma(n, m, s, tau):

	#print "sigma_fisher = ", 1./sqrt(fisher), "sig_asimov", float(s)/sqrt(2.*logLambda(n,m,s,tau))
	#from asimov approach
	if s>0 and logLambda(n,m,s,tau) > 0:
		#print "sigma from qmu"
		return float(s)/sqrt(2.*logLambda(n,m,s,tau))
	else:
		#print "sigma from fisher", s, logLambda(n,m,s,tau), 1/fisher(n,m,s,tau)
		return 1./fisher(m*tau+s,m,s,tau)


def F(qmu, mu, muprime, sigma):
	#print "debug mu=%f mu'=%f, sigma=%f, qmu=%f, dist_to_boundary = %e" %(mu, muprime, sigma, qmu, qmu-mu*mu/sigma/sigma )

	if qmu<mu*mu/sigma/sigma :
		#print "s = ", mu, "using first part", qmu, " ", sigma
		#print "debug ", sqrt(qmu), (mu-muprime)/sigma
		return norm.cdf(sqrt(qmu)-(mu-muprime)/sigma)
	else:
		#print "s = ", mu, "using second part", qmu, " ", sigma
		return norm.cdf( (qmu-(mu*mu-2.*mu*muprime)/sigma/sigma) / (2.*mu/sigma) )


def CLs_qmu(qmu, s, sig):
	#print "double scalar?", qmu, s, sig, (1.-F(qmu,s,s,sig)), (1.-F(qmu,s,0.0,sig))
	return (1.-F(qmu,s,s,sig))/(1.-F(qmu,s,0.0,sig))


def CLs(n, m, s, tau):
	#sig = sigma(n,m,s,tau)
	sig = sigma(n,m,s,tau)
	qmu = 2.*logLambda(n,m,s,tau)
	if shat(n,m,tau) > s :
		print "hit boundary"
		qmu=0

	#print "\ns = ", s, "qmu = ", qmu, "cutoff", s*s/sig/sig
	#print "sigma = ", sig
	#print "CLb = ", (1-F(qmu,s,0.,sig))
	#print "CLsb = ", F(qmu,s,s,sig)
	return (1.-F(qmu,s,s,sig))/(1.-F(qmu,s,0.0,sig))

def testReload():
	print "4"

def CLsArgumentWrapper(s,n,m,tau):
	return CLs(float(n),float(m),float(s),float(tau))-0.05

def CLsArgumentWrapper_qmu(s, bExp, deltaB, sigmaBand):
	tau = bExp/deltaB/deltaB
	sig = sigma(bExp,bExp*tau,s,tau)
	qmu = qmuBand(s, bExp, deltaB, sigmaBand)
	return CLs_qmu(float(qmu),float(s),float(sig))-0.05

def CLbArgumentWrapper(qmu, s,n,m,tau,sigmaBand):
	#using this to find value of q_mu that for sigmaBand of the s=0 distribution as a function of mu=s
	clbForBand = 1.-norm.cdf( sigmaBand)
	sig = sigma(n,m,s,tau)
	#sig = sigma(n,m,0.,tau)
	if s==0: #qmu is a delta function when mu=0
		print "qmu is a delta function when mu=0"
		return 1.;
	return F(float(qmu),float(s),float(0),float(sig))-clbForBand

def CLb(n,m,s,tau):
		#sig = sigma(n,m,s,tau)
	sig = sigma(n,m,s,tau)
	qmu = 2.*logLambda(n,m,s,tau)
	if shat(n,m,tau) >= s :
		print "hit boundary"
		qmu=0

	if s==0: #qmu is a delta function when mu=0
		print "qmu is a delta function when mu=0"
		return 1.;
	return F(float(qmu),float(s),float(0),float(sig))

def qmuBand(s, bExp, deltaB, sigmaBand) :
	# the value of q_mu that corresponds to the (-2,-1,0,1,2)sigmaBand of b-only for a hypothesized mu=s
	tau = bExp/deltaB/deltaB
	sig = sigma(bExp,bExp*tau,s,tau)
	xhi = (3*s/sig)**2 # a reasonable guess, but it breaks down for very small s near 0
	if xhi<0.01: xhi=0.01  # hack patch for small s
	xlo = xhi/100.
	try:
		qmuBand = brentq(CLbArgumentWrapper, xlo, xhi, args=(s,bExp,bExp*tau,tau,sigmaBand))
		#print "Approximate expected 95% CLs upper-limit = ", s 
		return qmuBand
	except:
		if F(float(0),float(s),float(0),float(sig)) > norm.cdf(sigmaBand):
			print "\n   brentq failed to find qmuBand, using a simple scan", xlo, xhi
			print "   P(delta)", F(float(0),float(s),float(0),float(sig)) , "CLb target =", norm.cdf(sigmaBand)
			print "   CLbArgumentWrapper", (s,bExp,bExp*tau,tau,sigmaBand)

		qmuBand=xlo
		#print "\n\ts, bExp, bexp*tau, sigmaBand", s,bExp,bExp*tau,tau,sigmaBand
		#print 'clb(0) = ', CLbArgumentWrapper(0, s, bExp,bExp*tau,tau, sigmaBand)
		#print 'clb(0.001) = ', CLbArgumentWrapper(0.001, s, bExp,bExp*tau,tau, sigmaBand)
		while CLbArgumentWrapper(qmuBand, s, bExp,bExp*tau,tau, sigmaBand) < 0.001 :
			qmuBand+=0.01
		#print "\n\tApproximate qmuBand = ", qmuBand
		return qmuBand


def ExpectedLimit(bExp, deltaB) :

	tau = bExp/deltaB/deltaB
	xhi = 3*sqrt(bExp)+3*deltaB #a reasonable guess of what is larger than the upper limit
	try:
		s = brentq(CLsArgumentWrapper, 0.001, xhi, args=(bExp,bExp*tau,tau))
		#print "Approximate expected 95% CLs upper-limit = ", s 
		return s
	except:
		print "\n\tbrentq failed (boundaries?) using a simple scan"
		s=0.001
		while CLs(bExp,bExp*tau,s,tau) < 0.95 :
			s+=0.01
		print "\n\tApproximate expected 95% CLs upper-limit = ", s 
		return s
	# use http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html#scipy.optimize.brentq


def ExpectedLimitBand(bExp, deltaB, sigmaBand) :
	xhi = 3*sqrt(bExp)+3*deltaB #a reasonable guess of what is larger than the upper limit
	try:
		s = brentq(CLsArgumentWrapper_qmu, 0.001, xhi, args=(bExp, deltaB, sigmaBand))
		#print "Approximate expected 95% CLs upper-limit = ", s 
		return s
	except:
		print "\n\tbrentq failed (boundaries?) using a simple scan"
		s=0.001
		while abs(CLsArgumentWrapper_qmu(s, bExp, deltaB, sigmaBand)) < 0.001 :
			s+=0.01
		print "\n\tCLsArgumentWrapper_qmu", (s, bExp, deltaB, sigmaBand)
		print "\n\tApproximate ", sigmaBand, " 95% CLs upper-limit = ", s 
		return s
 



def ObservedLimit(n, bExp, deltaB) :

	tau = bExp/deltaB/deltaB

	xhi = 3*sqrt(bExp)+3*deltaB #a reasonable guess of what is larger than the upper limit
	try:
		s = brentq(CLsArgumentWrapper, 0.001, xhi, args=(n,bExp*tau,tau))
		#print "Approximate expected 95% CLs upper-limit = ", s 
		return s
	except:
		print "brentq failed (boundaries?) using a simple scan"
		s=0.01
		while CLs(n,bExp*tau,s,tau) < 0.95 :
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

def ExpectedSignificance2(s, b, deltaB):
	tau = b/deltaB/deltaB
	q0 = 2.*logLambda(n=s+b,m=tau*b,s=0.,tau=tau)
	print q0, tau
	return sqrt(q0)

def ObservedSignificance(n, b, deltaB):
	tau = b/deltaB/deltaB
	q0 = 2.*logLambda(n=n,m=tau*b,s=0.,tau=tau)
	print q0, tau
	return sqrt(q0)


def ObsExpAndBands_Limits(n,bExp,deltaB):
	ol = ObservedLimit(n,bExp,deltaB)
	el = ExpectedLimit(bExp,deltaB)
	el_m2 = ExpectedLimitBand(bExp,deltaB,-2.)
	el_m1 = ExpectedLimitBand(bExp,deltaB,-1.)
	el_p1 = ExpectedLimitBand(bExp,deltaB,1.)
	el_p2 = ExpectedLimitBand(bExp,deltaB,2.)
	return (ol,  el_m2, el_m1, el, el_p1, el_p2)

if __name__ == '__main__':

	print 'qmuband(-1) = ', qmuBand(10, 100., 1., -1.)
	print 'qmuband(0) = ', qmuBand(10, 100., 1., 0.)
	print 'qmuband(+0.5) = ', qmuBand(10, 100., 1., 0.5)
	print 'qmuband(+0.8) = ', qmuBand(10, 100., 1., 0.8)
	print 'qmuband(+0.9) = ', qmuBand(10, 100., 1., 0.9)
	print 'qmuband(+0.95) = ', qmuBand(10, 100., 1., 0.95)
	print 'qmuband(1) = ', qmuBand(10, 100., 1., 1.)


	print '-2*logLambda = ', 2*logLambda(100, 100.*100., 10., 100)
	print '-2*logLambda(93,100) = ', 2*logLambda(93, 100.*100.-sqrt(100.*100./2.), 10., 100) #rough approx

	print 'expected limit ', ExpectedLimit(100, 1.) 
	print 'expected limit -1',ExpectedLimitBand(100, 1., -1.) 
	print 'expected limit 0',ExpectedLimitBand(100., 1., 0.) 
	print 'expected limit +0.5',ExpectedLimitBand(100, 1., 0.5) 
	#print 'expected limit +0.3',ExpectedLimitBand(100, 1., 0.3) 
	#print 'expected limit +0.5',ExpectedLimitBand(100, 1., 0.5) 
	#print 'expected limit +0.7',ExpectedLimitBand(100, 1., 0.7) 
	print 'expected limit +1',ExpectedLimitBand(100, 1., 1.) 

	#careful if shat(nsigma)> s

	print "running some tests"
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

	if	12.2 < sigma(100,50,50,1) and sigma(100,50,50,1) <12.3 :
		print  "ok" 
	else:
		print  "oops" 

	if	logLambda(100,50,50,1)==0 :
		print  "ok" 
	else:
		print  "oops" 

	if	ExpectedSignificance(50,100,10)>3. and ExpectedSignificance(50,100,10)<3.2 :
		print  "ok", "approximate significance for s=50,b=100+/-10 is",  ExpectedSignificance(50,100,10)
	else:
		print  "oops" 

	if	ExpectedLimit(50,7)>19.7 and ExpectedLimit(50,7)<19.8 :
		print  "ok",  "approximate limit for 50+/-7 is", ExpectedLimit(50,7)
	else:
		print  "oops" 

	if	ExpectedLimit(100,.1)>20. and ExpectedLimit(100,.1)<21 :
		print  "ok" , "approximate limit for 100+/-0.1 is", ExpectedLimit(100,0.1)
	else:
		print  "oops" 

	if	ExpectedLimit(50,50)>51 and ExpectedLimit(50,50)<53 :
		print  "ok", "approximate limit for 50+/-50 is", ExpectedLimit(50,50) 
	else:
		print  "oops" 

