// Author: Kyle Cranmer   October 2014
/*************************************************************************
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
  * you expect 50 +/- 7 background events
  * ExpectedLimit(50,7)
  * returns s_95 = 19.7 events

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
after making the replacements `n=bExp` and `m=bExp*tau`. 

 *************************************************************************/


#ifndef ROOT_Math_DistFuncMathCore
#include"Math/DistFuncMathCore.h"
#endif
#include <iostream>

double bhathat(double n, double m, double s, double tau){
	return (m+n-s-s*tau + sqrt(4*m*s*(1+tau) + pow(s+s*tau-m-n,2) ) ) / (2 * (1+tau )) ;
}

double shat(double n, double m, double tau){
	return (n*tau -m )/tau;
}

double bhat(double n, double m, double tau){
	return m/tau;
}


double logL(double n, double m, double s, double b, double tau){
	return (s+b) - n * log(s+b) + (tau*b) - m * log(tau*b);
}

double logLambda(double n, double m, double s, double tau){
	return logL(n, m, s, bhathat(n,m,s,tau), tau) 
		- logL(n,m,shat(n,m,tau), bhat(n,m,tau), tau);
}

double sigma(double n, double m, double s, double tau){
	double fisher = 8*m*n*pow(1+tau,2)* 
		(m*m+2*m*n+n*n+m*s-n*s+m*s*tau-n*s*tau+(n+m)*sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2))) 
		/sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)) 
		/pow(m+n-s-s*tau+sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)),2) 
		/pow(m+n+s+s*tau+sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)),2) ;
	//from asimov approach
	if(s>0 && logLambda(n,m,s,tau) > 0 )
		return s/sqrt(2.*logLambda(n,m,s,tau));
	else
		return 1./sqrt(fisher)	;

}

double F(double qmu, double mu, double muprime, double sigma){
	if(qmu < mu*mu/sigma/sigma){
		//print "s = ", mu, "using first part", qmu, " ", sigma
		//print "debug ", sqrt(qmu), (mu-muprime)/sigma
		return ROOT::Math::normal_cdf(sqrt(qmu)-(mu-muprime)/sigma);
	} else{
		//print "s = ", mu, "using second part", qmu, " ", sigma
		return ROOT::Math::normal_cdf( (qmu-(mu*mu-2*mu*muprime)/sigma/sigma) / (2*mu/sigma) );
	}
}

double CLs(double n, double m, double s, double tau){
	double qmu = 2*logLambda(n,m,s,tau);
	double sig = sigma(n,m,s,tau);
	if(shat(n,m,tau) > s)
		qmu=0.;

	return (1.-F(qmu,s,s,sig))/(1.-F(qmu,s,0,sig));
}

// function class with a member function for root finder
struct CLsHelper { 
	double n=0;
	double m=0;
	double tau=0;
	double evalCLs (double s ) const { 
   		return CLs(n,m,s,tau);
	}
};


double ExpectedLimit(double bExp, double deltaB) {

	double tau = bExp/deltaB/deltaB;

	double s=0.1;
	while(CLs(bExp,bExp*tau,s,tau) > 0.05 )
		s+=0.01;
	//cout << "Approximate expected 95% CLs upper-limit = " << s << endl;

	/*
	// to do, add a real root finder algorithm here
	CLsHelper helper;
	helper.n=bExp;
	helper.m=bExp*tau;
	helper.tau=tau;

	ROOTMathFunctor1D thisFunc(helper, &(helper.evalCLs));

 	ROOTMathFunctor1D f1D(&thisFunc);

	ROOTMathBrentRootFinder brf;
	brf.SetFunction(f1D);
	*/
	return s;

}


double ExpectedSignificance(double s, double b, double deltaB){
	// same idea, but this is expected significance (in sigma)
	return sqrt(2*( (s+b)*log( ((s+b)*(b+deltaB*deltaB)) / (b*b+(s+b)*deltaB*deltaB) )
				-( b*b/deltaB/deltaB)*log(1.+ (s*deltaB*deltaB)/(b*(b+deltaB*deltaB)) ) 
				  )
			);
}

void tests(){

	if(	shat(100,50,1)==50)
		cout << "ok" << endl;
	else
		cout << "oops" << endl;

	if(	bhat(100,50,1)==50)
		cout << "ok" << endl;
	else
		cout << "oops" << endl;

	if(	bhathat(100,50,50,1)==50)
		cout << "ok" << endl;
	else
		cout << "oops" << endl;

	if(	12.2 < sigma(100,50,50,1) && sigma(100,50,50,1) <12.3)
		cout << "ok" << endl;
	else
		cout << "oops" << endl;

	if(	logLambda(100,50,50,1)==0)
		cout << "ok" << endl;
	else
		cout << "oops" << endl;

	if(ExpectedSignificance(50,100,10)>3. && ExpectedSignificance(50,100,10)<3.2)
		cout <<   "ok" << "approximate significance for s=50,b=100+/-10 is" <<  ExpectedSignificance(50,100,10) << endl;
	else
		cout <<   "oops" << endl; 

	if(ExpectedLimit(50,7)>19.7 && ExpectedLimit(50,7)<19.8 )
		cout <<   "ok " <<  "approximate limit for 50+/-7 is " << ExpectedLimit(50,7) << endl;
	else
		cout <<   "oops" << endl; 

	if(	ExpectedLimit(100,.1)>20. && ExpectedLimit(100,.1)<21 )
		cout <<   "ok " << "approximate limit for 100+/-0.1 is " << ExpectedLimit(100,0.1) << endl;
	else
		cout <<   "oops" << endl; 

	if(	ExpectedLimit(50,50)>51 && ExpectedLimit(50,50)<53 )
		cout <<   "ok " << "approximate limit for 50+/-50 is " << ExpectedLimit(50,50)  << endl;
	else
		cout <<   "oops" << endl; 

}

