#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <iostream>
#include <string>
#include <vector>
#include <TH1D.h>
#include <TF1.h>
#include <TRandom3.h>

using namespace std;

class Experiment {

  public: 

    Experiment (double s, double b, double tau, TRandom3* ran);
    void generateData(double s, double b, double tau);
    double lnL(int n, int m, double s, double b, double tau);
    double lnL(int n, int m);
    double q0(int n, int m);
    double q0();
    double Z0();
    double s() { return m_s; }
    double b() { return m_b; }
    double tau() { return m_tau; }
    int n() { return m_n; }
    int m() { return m_m; }
    int Ntot(double sigrel);

  private:

    int       m_n;
    int       m_m;
    double    m_nd;
    double    m_md;
    double    m_s;
    double    m_b;
    double    m_tau;
    TRandom3* m_ran;

};

#endif
