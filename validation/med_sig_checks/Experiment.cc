#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <TMath.h>
#include <TF1.h>
#include <TRandom3.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/QuantFuncMathCore.h>
#include "Experiment.h"

using namespace std;

Experiment::Experiment(double s, double b, double tau, TRandom3* ran){
  m_s = s;
  m_b = b;
  m_tau = tau;
  m_n = -1;        // negative data value indicates not yet generated
  m_m = -1;
  m_nd = -1.;
  m_md = -1.;
  m_ran = ran;
}

void Experiment::generateData(double s, double b, double tau){
  m_n = m_ran->Poisson(s+b);
  m_nd = static_cast<double>(m_n);
  m_m = m_ran->Poisson(tau*b);
  m_md = static_cast<double>(m_m);
}

double Experiment::lnL(int n, int m, double s, double b, double tau){
  double nd = static_cast<double>(n);
  double md = static_cast<double>(m);
  double A = 0.;
  if ( n != 0 ) { A = nd*log(s + b); }
  double B = 0.;
  if ( m != 0 ) { B = md*log(tau*b); }
  double logL = A - (s + b) + B - tau*b;
  return logL;
}

double Experiment::lnL(int n, int m){   // uses object's values of s, b, tau
  double s = m_s;
  double b = m_b;
  double tau = m_tau;
  return this->lnL(n, m, s, b, tau);
}

double Experiment::q0(int n, int m){
  double s = m_s;
  double b = m_b;
  double tau = m_tau;
  double nd = static_cast<double>(n);
  double md = static_cast<double>(m);
  double bHatHat0 = (nd + md)/(1. + m_tau);
  double bHat = md/m_tau;  
  double sHat = nd - bHat;
  double lnLambda0 = this->lnL(n, m, 0, bHatHat0, tau) - 
                     this->lnL(n, m, sHat, bHat, tau);
  double qZero = 0.;
  if ( sHat > 0. ) {
    qZero = - 2.*lnLambda0;
  }
  return qZero;
}

double Experiment::q0(){
  int n = m_n;
  int m = m_m;
  return this->q0(n, m);
}

int Experiment::Ntot(double sigrel){

  // sigrel = specified relative accuracy for Z
  // return value is number of MC iterations needed to achieve this

  const double pi = 3.14159265;
  const double rootTwoPi = sqrt(2.*pi);
  double q0Obs = this->q0();    // based on m_n and m_m
  double Z = sqrt(q0Obs);       // Asymptotic value
  double p = 1. - TMath::Freq(Z);
  double u = -2.*log(rootTwoPi*p);
  double dZdp = (1 - u)/(u*p*Z);
  double N = dZdp*dZdp*p*(1.-p)/(Z*Z*sigrel*sigrel);
  int ntot = static_cast<int>(N);
  // cout << "q0obs, ntot = " << q0Obs << "  " << ntot << endl;
  return ntot;
}

double Experiment::Z0(){

  const double pi = 3.14159265;
  const double rootTwoPi = sqrt(2.*pi);

  // computes p-value of s=0 using objects vales
  // of m_n and m_m

  int ntotMin = 1000;
  int nq0HigherMin = 5;
  int ntotMax = 40000000;
  double ZCut = 5.5;       // for ZAsymp > ZCut, use ZAsymp

  double q0Obs = this->q0();    // based on m_n and m_m
  double ZAsymp = sqrt(q0Obs);
  double Z = 999.;

  // cout << "n, m, ZAsymp" << "  " << m_n << "  " << m_m << "  " 
  //      << ZAsymp << endl;

  if ( m_n == 0 ) {
    Z = 0.;             // no discovery if n = 0
  }
  else if ( ZAsymp >= ZCut ) {
    Z = ZAsymp;         //  this should not happen more than half the time
  }
  else {
    int nq0Higher = 0;
    // bool stopLoop = false;

    // estimate how many iterations needed to find Z to within
    // a specified relative uncertainty, sigrel.

    double sigrel = 0.01;
    int ntot = this->Ntot(sigrel);
    // cout << "n, m, ZAsymp, ntot = " << m_n << "  " << m_m << "  " << 
    //  ZAsymp << "  " << ntot << endl;

    if ( ntot < ntotMin ) { ntot = ntotMin; }
    if ( ntot > ntotMax ) { ntot = ntotMax; }

    //    while ( !stopLoop ) {
    for (int i=0; i<ntot; i++){
      int n = m_ran->Poisson(m_b);       // bkg-only
      int m = m_ran->Poisson(m_tau*m_b);
      double qZero = this->q0(n, m);
      if ( qZero >= q0Obs ) {  nq0Higher++; }


      /// ntot++;
      // double sigmaZ_over_Z = 999.;
      // if ( ntot%1000 == 0 && nq0Higher >= 1 ) {
      //   double p = static_cast<double>(nq0Higher)/static_cast<double>(ntot);
      //   double u = -2.*log(rootTwoPi*p);
      //   double Z = sqrt(u - log(u));
      //   double Vp =  p*(1.-p) / static_cast<double>(ntot);
      //   double dZdp = (1 - u)/(u*p*Z);
      //   double VZ = dZdp*dZdp*Vp;
      //   sigmaZ_over_Z = sqrt(VZ)/Z;

    }
    //      stopLoop = sigmaZ_over_Z < 0.02;

      // stopLoop = (ntot >= ntotMin && nq0Higher >= nq0HigherMin)
      //		   || ntot >= ntotMax;

    double one_minus_p0 = static_cast<double>(ntot-nq0Higher) / 
              static_cast<double>(ntot);

    Z = 0.;
    if ( one_minus_p0 >= 1. ) {
      Z = 999.;
    }
    else if ( one_minus_p0 <= 0. ) {
      Z = -999.;
    }
    else {
      Z = ROOT::Math::normal_quantile(one_minus_p0, 1.);
    }

  }

// cout << "ntot, nq0Higher = " << ntot << "  " << nq0Higher << endl;
  // cout << "returning Z = " << Z << endl;

  return Z;

}
