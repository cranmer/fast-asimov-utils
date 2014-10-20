// Author: Kyle Cranmer   October 2014
/*************************************************************************
 * Copyright (C) Kyle Cranmer
_ License LGPL v2.1 (for ROOT compatibility, happy to make it BSD for other purposes) _

This code quickly calculates the expected 95% CLs upper limit on the number of
signal events, s, given an expected background, bExp, and 
uncertainty on the background estimate, deltaB.
The background uncertainty is absolute (not relative) and is uncertainty
on the mean background (so you don't include Poisson fluctuatiosn in this number).
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

double sigma(double n, double m, double s, double tau){
	return 8*m*n*pow(1+tau,2)*
	(m*m+2*m*n+n*n+m*s-n*s+m*s*tau-n*s*tau+(n+m)*sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)))
	/sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2))
	/pow(m+n-s-s*tau+sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)),2)
	/pow(m+n+s+s*tau+sqrt(4*m*s*(1+tau)+pow(m+n-s*(1+tau),2)),2);	

}

double logL(double n, double m, double s, double b, double tau){
	return (s+b) - n * log(s+b) + (tau*b) - m * log(tau*b);
}

double logLambda(double n, double m, double s, double tau){
	return logL(n, m, s, bhathat(n,m,s,tau), tau) 
		- logL(n,m,shat(n,m,tau), bhat(n,m,tau), tau);
}

double F(double qmu, double mu, double muprime, double sigma){
	return ROOT::Math::normal_cdf(sqrt(qmu)-(mu-muprime)/sigma);
}

double CLs(double n, double m, double s, double tau){
	double qmu = 2*logLambda(n,m,s,tau);
	double sig = sigma(n,m,s,tau);
	if(shat(n,m,tau) > s)
		return 0.5;

	return F(qmu,s,s,sig)/(1-F(qmu,s,0,sig));
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
	while(CLs(bExp,bExp*tau,s,tau) < 0.95 )
		s+=0.01;
	cout << "Approximate expected 95% CLs upper-limit = " << s << endl;

/*
	CLsHelper helper;
	helper.n=bExp;
	helper.m=bExp*tau;
	helper.tau=tau;

	ROOT::Math::Functor1D thisFunc(helper, &(helper.evalCLs));

 	ROOT::Math::Functor1D f1D(&thisFunc);

	ROOT::Math::BrentRootFinder brf;
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

	if(	149 < 1/sigma(100,50,50,1) && 1/sigma(100,50,50,1) <151)
		cout << "ok" << endl;
	else
		cout << "oops" << endl;

	if(	logLambda(100,50,50,1)==0)
		cout << "ok" << endl;
	else
		cout << "oops" << endl;

}

