// Program to illustrate simple statistical test
// Glen Cowan, RHUL Physics, June 2010

#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TMath.h>
#include "Experiment.h"

using namespace std;

int main(){

// Open output file (apparently needs to be done before booking)

  TFile* file = new TFile("medsigMC.root", "recreate");

// Book histograms

  TH1D* h_q0 = new TH1D("h_q0", "q0", 400, 0., 40.);

// Create a TRandom3 object to generate random numbers

  int seed = 12345;
  TRandom3* ran = new TRandom3(seed);

  int numSigval = 3;
  vector<double> sigVec(numSigval);
  sigVec[0] = 0.5;
  sigVec[1] = 1.0;
  sigVec[2] = 2.0;

  vector<double> relSigVec(numSigval);
  relSigVec[0] = 0.2;
  relSigVec[1] = 0.5;
  relSigVec[2] = 1.0;

  vector<double> tauVec(numSigval);
  tauVec[0] = 0.5;
  tauVec[1] = 1.0;
  tauVec[2] = 2.0;

  double s = 5.;

  double bMin = 0.1;
  double bMax = 100;
  double logbMin = log(bMin)/log(10);
  double logbMax = log(bMax)/log(10);
  int numBval = 20;

  vector<double> Z_naive(numSigval);
  vector<double> Z_wald(numSigval);
  vector<double> meanZ0(numSigval);         // sqrt(q0)
  vector<double> medZ0(numSigval);          // sqrt(q0)
  vector<double> medZ0MC(numSigval);        // exact
  vector<double> Z_bi(numSigval);           // Asimov w/ binomial

  for (int i=0; i<numBval; i++){
    double logb = (logbMax - logbMin)*static_cast<double>(i) / 
                  static_cast<double>(numBval-1)  +  logbMin;
    double b = pow(10., logb);

    for (int j=0; j<numSigval; j++){

      //       double sig = sigVec[j];

      double sig = relSigVec[j] * b;
      double tau = 0;
      if ( sig > 0 ) { tau = b / (sig*sig); }

      // double tau = tauVec[j];
      // double sig = sqrt(b/tau);

      Z_naive[j] = s/sqrt(b + sig*sig);
      double p_bi = TMath::BetaIncomplete(1./(1.+tau), s+b, tau*b + 1.);
      Z_bi[j] = TMath::NormQuantile(1. - p_bi);
      if ( sig == 0 ) {
        Z_wald[j] = sqrt( 2.*((s + b)*log(1. + s/b) - s) );
      }
      else {
        double X = (s+b)*log( (s+b)*(b+sig*sig)/(b*b + (s+b)*sig*sig) );
        double Y = (b*b/(sig*sig))*log(1. + sig*sig*s/(b*(b+sig*sig)));
        Z_wald[j] = sqrt(2.*(X - Y));
      }

      Experiment exper(s, b, tau, ran);
      int numExp = 101;
      int middle = (numExp - 1)/2;
      vector<double> Z0Vec;
      vector<double> Z0MCVec;
      for (int k=0; k<numExp; ++k){
        exper.generateData(s, b, tau);     
        double q0 = exper.q0();       
        double Z0 = 0.;
        if ( q0 > 0. ) { Z0 = sqrt(q0); }

        double Z0MC = exper.Z0();
	// cout << s << "  " << b << "  " << exper.n() << "  " << exper.m() 
        //     << "  " << Z0 << "  " << Z0MC << endl;


        if ( j == 3 && i == 40 ) {
          h_q0->Fill(q0);
          int n = exper.n();
          int m = exper.m();
          double bHat = static_cast<double>(m)/tau;
	  // cout << b << "  " << tau << "  " << bHat << "  "
          //     << n << "  " << m << "  " << q0 << "  " << Z0 << endl;
        }
        Z0Vec.push_back(Z0);
        Z0MCVec.push_back(Z0MC);
        // if ( k%10 == 0 ) { cout << "generated exp " << k+1 << endl; }
      }
      sort(Z0Vec.begin(), Z0Vec.end());
      medZ0[j] = Z0Vec[middle];
      sort(Z0MCVec.begin(), Z0MCVec.end());
      medZ0MC[j] = Z0MCVec[middle];




    }

    cout << b << "  " << Z_naive[0] << "  " <<
    Z_naive[1] << "  " << Z_naive[2] << "  " << 
    Z_wald[0] << "  " << Z_wald[1] << "  " << 
    Z_wald[2] << "  " << medZ0[0] << "  " << 
    medZ0[1] << "  " << medZ0[2] << "  " << 
    medZ0MC[0] << "  " << medZ0MC[1] << "  " << 
    medZ0MC[2] << "  " << Z_bi[0] << "  " <<
    Z_bi[1] << "  " << Z_bi[2] << "  " << endl;

  }

// Store all histograms in the output file and close up

  file->Write();
  file->Close();

  return 0;
}
