/* 
 *  FrequencyEstimatorTest.cpp
 *  FrequencyEstimatorTest
 *
 *    This is a simple test function for the class FrequencyEstimator.
 *
 *  Created by Andre Pollok on 08/02/2013.
 *  Copyright (c) 2013 ITR, University of South Australia. All rights reserved.
 */

#include <cstdlib>
#include <iostream>
#include <complex>
#include <vector>
#include "FrequencyEstimator.h"

using namespace std;

// tolerance for tests
VALTYPE tol = 1e-5;

VALTYPE Ts = 1.0;   // sample period

vector<complexT> txSamp; // vector of transmitted samples
vector<complexT> rxSamp; // vector of received samples
vector<VALTYPE>  tSamp;  // vector of sample times
VALTYPE tauref = 10*Ts; // time reference point
unsigned int numSamp = 100;
unsigned int numRepeat = 1; // number of times the estimator is run (numRepat>1 only useful for profiling)

// channel parameter offsets
VALTYPE   ftrue      = -0.813/(2*Ts); // frequency offset in Hz
VALTYPE   frtrue     = +0.0009; // frequency rate offset in Hz/s
complexT  cgaintrue  = (complexT)std::polar(0.5, M_PI_4);

// estimator and parameter estimates
FrequencyEstimator* est;
unsigned int Nfft = 128; // FFT size
VALTYPE frmin = -0.001; // lower limit for frequency rate search range in Hz/s
VALTYPE frmax = +0.001; // upper limit for frequency rate search range in Hz/s
VALTYPE search_oversample = 1; // oversampling factor of the frequency rate search
VALTYPE   fhat;     // frequency estimate
VALTYPE   frhat;    // frequency rate estimate
complexT  cgainhat; // complex gain estimate

int main(int argc, char** argv) {

  for (int m = 0; m < numSamp; m++) tSamp.push_back(m*Ts-tauref); // vector of sample times
  for (int m = 0; m < numSamp; m++) txSamp.push_back(1); // assume Tx signal is a unit-amplitude tone at DC
  
  cout << "Testing estimator with frequency offset..." << "\t";
  cout.flush();
  for (int m = 0; m < numSamp; m++) {
    rxSamp.push_back(txSamp[m] * cgaintrue * (complexT)std::polar<VALTYPE>(1.0, 2.0*M_PI*ftrue*tSamp[m]));
  }
  // construct and run estimator
  est = new FrequencyEstimatorLeastSquares(Nfft,Ts,tauref);
  for (unsigned int i=0; i<numRepeat; i++) {
//    est->estimate(rxSamp, txSamp);
    est->estimate(rxSamp);
    est->refineEstimate();
  }
  fhat      = est->frequencyEstimate();
  cgainhat  = est->complexGainEstimate();
  //cout << fhat - ftrue << cgainhat - cgaintrue << endl;
  bool pass = true;
  if (abs(ftrue     - fhat)     > tol) pass = false;
  if (abs(cgaintrue - cgainhat) > tol) pass = false;
  if (pass) std::cout << "PASS" << std::endl;
  else std::cout << "FAIL" << std::endl;
  
  cout << "Testing estimator with frequency and frequency rate offsets..." << "\t";
  cout.flush();
  for (int m = 0; m < numSamp; m++) {
    rxSamp[m] = txSamp[m] * cgaintrue * (complexT)std::polar<VALTYPE>(1.0, 2.0*M_PI * (ftrue*tSamp[m] + (0.5*frtrue)*tSamp[m]*tSamp[m]));
  }
  // construct and run estimator
  est = new FrequencyEstimatorLeastSquaresWithRate(Nfft,Ts,tauref,frmin,frmax,search_oversample);
  for (unsigned int i=0; i<numRepeat; i++) {
//    est->estimate(rxSamp, txSamp);
    est->estimate(rxSamp);
    est->refineEstimate();
  }
  fhat      = est->frequencyEstimate();
  frhat     = est->frequencyRateEstimate();
  cgainhat  = est->complexGainEstimate();
  pass = true;
  if (abs(ftrue     - fhat)     > tol) pass = false;
  if (abs(frtrue    - frhat)    > tol) pass = false;
  if (abs(cgaintrue - cgainhat) > tol) pass = false;
  if (pass) std::cout << "PASS" << std::endl;
  else std::cout << "FAIL" << std::endl;
  
  delete est;
  
  return 0;
}
