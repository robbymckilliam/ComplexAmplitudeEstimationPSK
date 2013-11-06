/* 
 * File:   OerderMeyerMasseyTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 20/02/2013, 4:55:35 PM
 */

#include "mex.h"
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <complex>
#include "FinitePulse.h"
#include "Util.h"
#include "OerderMeyerMassey.h"
#include "FastAlgorithm.h"

#ifndef VALTYPE
#define VALTYPE VALTYPE
#endif
#ifndef MSGTYPE
#define MSGTYPE VALTYPE
#endif

using namespace std;
using namespace TimeOffset;


static VALTYPE tol = 1e-2;
static VALTYPE T = 1.0; //symbol period
static unsigned int p = 3;
static unsigned int q = 17;
static unsigned int c = 5; //new c, its will divide T
static VALTYPE Ts = p * T / q; //sample period
static VALTYPE Delta = T/c;
static VALTYPE taumin = 10.0;
static VALTYPE taumax = 30.0;
static unsigned int L = 100;
static unsigned int N = (int) ceil((T * L + taumax) / Ts);
static int numzeros = 2;
static vector<int> P;
static vector<int> D;
static vector<complex<VALTYPE> > pilots;
static vector<complex<VALTYPE> > s;
static vector<complex<VALTYPE> > r;
static unsigned int M = 4; //QPSK
static int numpilots = 30;

    //transmitted signal
complex<VALTYPE> x(vector<complex<VALTYPE> >&s, TruncatedRootRaisedCosine& tx, VALTYPE T, int L, VALTYPE t) {
      int mini = max(1, (int)ceil((t - tx.tmax())/T));
      int maxi = min(L, (int)floor((t - tx.tmin())/T));
      complex<VALTYPE> sum(0,0);
      for(int i = mini; i <= maxi; i++) sum += s[i-1] * tx.pulse(t - i*T);
      return sum;
};

void testwithstore() {
  //generator.seed(10);
  const VALTYPE benchtime = 10.0; //number of seconds to run each benchmark for.
  const VALTYPE tol = 1e-2;

  const unsigned int M = 4;
  const unsigned int c = 4;
  const VALTYPE T = 1.0; //symbol period
  const unsigned int p = 3;
  const unsigned int q = 17;
  const VALTYPE Ts = p*T/q; //sample period
  const VALTYPE taumin = 11.5;
  const VALTYPE taumax = 110.0;
  
  const VALTYPE tau0 = (taumax - taumin)/2 + taumin;
  const complex<VALTYPE> a0 = polar<VALTYPE>(1,0.1*pi);
  const int L = 100;
 
   //const unsigned int numpilots = L/10; //about 10% pilots 
    const unsigned int numpilots = 30;
    //symbol position setup
    vector<int> P;
    for(int i = 0; i < numpilots; i++) P.push_back(i); //pilots at the front
    vector<int> D;
    for(int i = numpilots; i < L; i++) D.push_back(i); //data at the back
    vector<complex<VALTYPE> > s;
    for (int m = 1; m <= L; m++) s.push_back(polar<VALTYPE>(1.0, 2 * pi * (rand()%M) / M));
    vector<complex<VALTYPE> > pilots;
    for (int m = 0; m < numpilots; m++) pilots.push_back(s[m]);

    TruncatedRootRaisedCosine tx(T,1.0/2.0,10);
    OerderMeyerAndMassey<TruncatedRootRaisedCosine> est(P,D,pilots,tx,T,Ts,taumin,taumax);
    //Direct<TruncatedRootRaisedCosine> est(P,D,pilots,tx,T,Ts,taumin,taumax,c,p,q);
  
        //number of samples
    unsigned int N = (unsigned int) ceil((T*(L+40)+taumax)/Ts); //number of samples
    
        //sampled received signal
    vector<complex<VALTYPE> > r;
    for(int n = 1; n <= N; n++) r.push_back(a0*x(s,tx,T,L,n*Ts-tau0));
    
    VALTYPE tauhat = est.estimate(r);
    
    std::cout << "OerderMeyerAndMassey test ... ";
    bool pass = fabs(tau0 - tauhat) < tol;
    if (pass) std::cout << "PASS " << std::endl;
    else std::cout << "FAIL " << tau0 << ", " << tauhat << std::endl;
    
}

void testtabluated() {
  //generator.seed(10);
  const VALTYPE benchtime = 10.0; //number of seconds to run each benchmark for.
  const VALTYPE tol = 1e-2;

  const unsigned int M = 4;
  const unsigned int c = 4;
  const VALTYPE T = 1.0; //symbol period
  const unsigned int p = 3;
  const unsigned int q = 17;
  const VALTYPE Ts = p*T/q; //sample period
  const VALTYPE taumin = 11.5;
  const VALTYPE taumax = 110.0;
  
  const VALTYPE tau0 = (taumax - taumin)/2 + taumin;
  const complex<VALTYPE> a0 = polar<VALTYPE>(1,0.1*pi);
  const int L = 100;
  
   //const unsigned int numpilots = L/10; //about 10% pilots 
    const unsigned int numpilots = 30;
    //symbol position setup
    vector<int> P;
    for(int i = 0; i < numpilots; i++) P.push_back(i); //pilots at the front
    vector<int> D;
    for(int i = numpilots; i < L; i++) D.push_back(i); //data at the back
    vector<complex<VALTYPE> > s;
    for (int m = 1; m <= L; m++) s.push_back(polar<VALTYPE>(1.0, 2 * pi * (rand()%M)/M) );
    vector<complex<VALTYPE> > pilots;
    for (int m = 0; m < numpilots; m++) pilots.push_back(s[m]);

    TruncatedRootRaisedCosine tx(T,1.0/2.0,10);
    TabulatedOerderMeyerAndMassey<TruncatedRootRaisedCosine> est(P,D,pilots,tx,T,Ts,taumin,taumax);
    //Direct<TruncatedRootRaisedCosine> est(P,D,pilots,tx,T,Ts,taumin,taumax,c,p,q);
  
        //number of samples
    unsigned int N = (unsigned int) ceil((T*(L+40)+taumax)/Ts); //number of samples
    
        //sampled received signal
    vector<complex<VALTYPE> > r;
    for(int n = 1; n <= N; n++) r.push_back(a0*x(s,tx,T,L,n*Ts-tau0));
    
    VALTYPE tauhat = est.estimate(r);
    
    std::cout << "Tabulated OerderMeyerAndMassey test ... ";
    bool pass = fabs(tau0 - tauhat) < tol;
    if (pass) std::cout << "PASS " << std::endl;
    else std::cout << "FAIL " << tau0 << ", " << tauhat << std::endl;
    
}

void testalignedtabluated() {
  //generator.seed(10);
  const VALTYPE benchtime = 10.0; //number of seconds to run each benchmark for.
  const VALTYPE tol = 1e-2;

  const unsigned int M = 4;
  const unsigned int c = 4;
  const VALTYPE T = 1.0; //symbol period
  const unsigned int p = 3;
  const unsigned int q = 17;
  const VALTYPE Ts = p*T/q; //sample period
  const VALTYPE taumin = 11.5;
  const VALTYPE taumax = 110.0;
  
  const VALTYPE tau0 = (taumax - taumin)/2 + taumin;
  const complex<VALTYPE> a0 = polar<VALTYPE>(1,0.1*pi);
  const int L = 100;
 
   //const unsigned int numpilots = L/10; //about 10% pilots 
    const unsigned int numpilots = 30;
    //symbol position setup
    vector<int> P;
    for(int i = 0; i < numpilots; i++) P.push_back(i); //pilots at the front
    vector<int> D;
    for(int i = numpilots; i < L; i++) D.push_back(i); //data at the back
    vector<complex<VALTYPE> > s;
    for (int m = 1; m <= L; m++) s.push_back(polar<VALTYPE>(1.0, 2 * pi * (rand()%M)/M) );
    vector<complex<VALTYPE> > pilots;
    for (int m = 0; m < numpilots; m++) pilots.push_back(s[m]);

    TruncatedRootRaisedCosine tx(T,1.0/2.0,10);
    AlignedTabulatedOerderMeyerAndMassey<TruncatedRootRaisedCosine> est(P,D,pilots,tx,T,Ts,taumin,taumax);
    //Direct<TruncatedRootRaisedCosine> est(P,D,pilots,tx,T,Ts,taumin,taumax,c,p,q);
  
        //number of samples
    unsigned int N = (unsigned int) ceil((T*(L+40)+taumax)/Ts); //number of samples

        //sampled received signal
    vector<complex<VALTYPE> > r;
    for(int n = 1; n <= N; n++) r.push_back(a0*x(s,tx,T,L,n*Ts-tau0));
    
    VALTYPE tauhat = est.estimate(r);
    
    std::cout << "Aligned Tabulated OerderMeyerAndMassey test ... ";
    bool pass = fabs(tau0 - tauhat) < tol;
    if (pass) std::cout << "PASS " << std::endl;
    else std::cout << "FAIL " << tau0 << ", " << tauhat << std::endl;
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    testwithstore();
    testtabluated();
    testalignedtabluated();
 
}

