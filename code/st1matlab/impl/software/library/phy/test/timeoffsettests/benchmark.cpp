/* 
 * File:   montecarlo.cpp
 * Author: Robby McKilliam
 */

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "Util.h"
#include "FinitePulse.h"
#include "TimeOffsetEstimator.h"
#include "Recursive.h"
#include "Direct.h"
#include "OerderMeyerMassey.h"
#include <random>
#include <ctime>

using namespace std;
using namespace TimeOffset;

int main(int argc, char** argv) {

  default_random_engine generator; 
  const double benchtime = 20.0; //number of seconds to run each benchmark for.

  const double iters = 4000;
  const unsigned int M = 4;
  const unsigned int c = 4;
  const double T = 1.0; //symbol period
  const unsigned int p = 1;
  const unsigned int q = 4;
  const double Ts = p*T/q; //sample period
  const double taumin = 40.0;
  const double taumax = 50.0;
  const double tau0 = 46.1;
  const complex<double> a0 = polar<double>(1,0.1*pi);
  //const vector<int> Ls = {11,15,20,30,40,50,60,80,100,120,140,160,180,200,250,300,400,500,600,800,1000,1200,1500,2000,2600,3500};
  const vector<int> Ls = {3990};
 

  //SNRs
  double SNRdB = 10.0;
  double SNR = pow(10.0, SNRdB/10.0);
  normal_distribution<double> noise(0.0,T/Ts/SNR/2);
  uniform_int_distribution<int> unifM(1, M); //for generating M-PSK symbols

  for(auto L : Ls) {

    const unsigned int numpilots = L/10; //about 10% pilots 
    //symbol position setup
    vector<int> P;
    for(int i = 0; i < numpilots; i++) P.push_back(i); //pilots at the front
    vector<int> D;
    for(int i = numpilots; i < L; i++) D.push_back(i); //data at the back
    vector<complex<double>> s;
    for (int m = 1; m <= L; m++) s.push_back(polar<double>(1.0, 2 * pi * unifM(generator) / M));
    vector<complex<double>> pilots;
    for (int m = 0; m < numpilots; m++) pilots.push_back(s[m]);

    //the transmit pulse
    //TruncatedRootRaisedCosine tx(T,1.0/2.0,4);
    FinitePulseWithLookupTable<TruncatedRootRaisedCosine> tx(TruncatedRootRaisedCosine(T,1.0/3.0,10),5000);
    //list of estimators we will run
    
    Recursive<FinitePulseWithLookupTable<TruncatedRootRaisedCosine>> rec(P,D,pilots,tx,T,Ts,taumin,taumax,c,p,q,1e-4);
    //Direct dir(P,D,pilots,new FinitePulseWithLookupTable(new TruncatedRootRaisedCosine(T,1.0/3.0,15),10000),T,Ts,taumin,taumax,c,p,q);
    OerderMeyerAndMassey<FinitePulseWithLookupTable<TruncatedRootRaisedCosine>> om(P,D,pilots,tx,T,Ts,taumin,taumax);
    const vector<TimeOffsetEstimator*> ests = {&om};

    //number of samples
    unsigned int N = (unsigned int) ceil((T*(L+40)+taumax)/Ts); //number of samples
    
    //transmitted signal
    auto x = [&s,&tx,T,L] (double t) {
      int mini = max(1, (int)ceil((t - tx.tmax())/T));
      int maxi = min(L, (int)floor((t - tx.tmin())/T));
      complex<double> sum(0,0);
      for(int i = mini; i <= maxi; i++) sum += s[i-1] * tx.pulse(t - i*T);
      return sum;
    };

      //sampled received signal
      vector<complex<double>> r;
	for(int n = 1; n <= N; n++) r.push_back(a0*x(n*Ts-tau0) + complex<double>(noise(generator), noise(generator)));
      

    for( auto est : ests ) {
      cout << "Benchmarking " << est->name() << " L = " << L << " ... ";
      long iters = 0;
      clock_t started = clock();
      while( ((double)(clock() - started))/CLOCKS_PER_SEC < benchtime) {
	est->estimate(r);
	iters++;
      }
      clock_t stopped = clock();
      cout << " requires  " << ((double)(stopped - started))/CLOCKS_PER_SEC/iters/L*1000 << " milliseconds per symbol" << endl;
    }

  }

}
