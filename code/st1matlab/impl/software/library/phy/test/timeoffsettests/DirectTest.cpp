/* 
 * File:   DirectTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 25 January 2013, 7:31 PM
 */

#include "mex.h"
#include <cstdlib>
#include <iostream>
#include <complex>
#include "FastAlgorithm.h"
#include "TabulatedFastAlgorithm.h"
#include "FinitePulse.h"

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

void testDirect() {

    Direct<TruncatedSincPulse> dir(P, D, pilots, TruncatedSincPulse(T, numzeros), T, Ts, taumin, taumax, c, p, q);
    Naive<TruncatedSincPulse> nai(P, D, pilots, TruncatedSincPulse(T, numzeros), T, Ts, taumin, taumax, c);
    VALTYPE tauhatdir = dir.estimate(r);
    VALTYPE tauhatnai = nai.estimate(r);
    
    std::cout << "Direct Zgrid test ... ";
    bool pass = true;
    vector<VALTYPE> zgrid = dir.getZgrid();
    for(int i = 0; i < zgrid.size(); i++ ) {
        //std::cout << zgrid[i] << ", " << nai.Z(taumin + i*Delta) << std::endl;
        pass &= fabs(zgrid[i] - nai.Z(taumin + i*Delta)) < tol;
        pass &= fabs(zgrid[i] - dir.Z(taumin + i*Delta)) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "Direct Ygrid test ... ";
    pass = true;
    vector<VALTYPE> ygrid = dir.getYgrid();
    for(int i = 0; i < ygrid.size(); i++ ) {
        pass &= fabs(ygrid[i] - std::abs(nai.Y(taumin + i*Delta))) < tol;
        pass &= fabs(ygrid[i] - std::abs(dir.Y(taumin + i*Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "Direct SSgrid test ... ";
    pass = true;
    vector<VALTYPE> ssgrid = dir.getSSgrid();
    for(int i = 0; i < ssgrid.size(); i++ ) {
        pass &= fabs(ssgrid[i] - std::abs(nai.SS(taumin + i*Delta))) < tol;
        pass &= fabs(ssgrid[i] - std::abs(dir.SS(taumin + i*Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "Direct coarse estimate... ";
    //std::cout << dir.coarseMaximiseSS() << ", " << nai.coarseMaximiseSS() << std::endl;
    pass = fabs(dir.coarseMaximiseSS() - nai.coarseMaximiseSS()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "Direct fine estimate... ";
    pass = fabs(tauhatdir - tauhatnai) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
}

void testTabulatedDirect() {
    
  TabulatedDirect<TruncatedSincPulse> dir(P, D, pilots, TruncatedSincPulse(T, numzeros), T, Ts, taumin, taumax, c, p, q, 10000);
    Naive<TruncatedSincPulse> nai(P, D, pilots, TruncatedSincPulse(T, numzeros), T, Ts, taumin, taumax, c);
    VALTYPE tauhatdir = dir.estimate(r);
    VALTYPE tauhatnai = nai.estimate(r);
    
    std::cout << "TabulatedDirect Zgrid test ... ";
    bool pass = true;
    vector<VALTYPE> zgrid = dir.getZgrid();
    for(int i = 0; i < zgrid.size(); i++ ) {
      //std::cout << zgrid[i] << ", " << nai.Z(taumin + i*Delta) << std::endl;
        pass &= fabs(zgrid[i] - nai.Z(taumin + i*Delta)) < tol;
        pass &= fabs(zgrid[i] - dir.Z(taumin + i*Delta)) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "TabulatedDirect Ygrid test ... ";
    pass = true;
    vector<VALTYPE> ygrid = dir.getYgrid();
    for(int i = 0; i < ygrid.size(); i++ ) {
        pass &= fabs(ygrid[i] - std::abs(nai.Y(taumin + i*Delta))) < tol;
        pass &= fabs(ygrid[i] - std::abs(dir.Y(taumin + i*Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "TabulatedDirect SSgrid test ... ";
    pass = true;
    vector<VALTYPE> ssgrid = dir.getSSgrid();
    for(int i = 0; i < ssgrid.size(); i++ ) {
        //std::cout << taumin + i*Delta << ", " << ssgrid[i] << ", " << std::abs(nai.SS(taumin + i*Delta)) << std::endl;
        pass &= fabs(ssgrid[i] - std::abs(nai.SS(taumin + i*Delta))) < tol;
        pass &= fabs(ssgrid[i] - std::abs(dir.SS(taumin + i*Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "TabulatedDirect coarse estimate... ";
    //std::cout << dir.coarseMaximiseSS() << ", " << nai.coarseMaximiseSS() << std::endl;
    pass = fabs(dir.coarseMaximiseSS() - nai.coarseMaximiseSS()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "TabulatedDirect fine estimate... ";
    pass = fabs(tauhatdir - tauhatnai) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    for (int m = 0; m < numpilots; m++) P.push_back(m);
    for (int m = numpilots; m < L; m++) D.push_back(m);
    for (int m = 1; m <= L; m++) s.push_back(polar<VALTYPE>(1.0, (2*pi*(rand()%M))/M));
    for ( int i = 0; i < P.size(); i++ ) pilots.push_back(s[P[i]]);
    for (int n = 1; n <= N; n++) r.push_back(polar<VALTYPE>(1.0, -0.1 * n));

    testDirect();
    testTabulatedDirect();
 
}

