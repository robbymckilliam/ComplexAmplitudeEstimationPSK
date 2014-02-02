/* 
 * File:   turbosyncsim.cpp
 * Author: Robby McKilliam
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "Util.h"
#include "CodedConstellation.h"
#include "CoherentMackenthun.h"
#include "LDPCDec.h"
#include "TurboSync.h"
#include <random>
#include <functional>
#include <future>
#include <chrono>

using namespace std;

typedef std::chrono::duration<double, std::ratio<1>> seconds;

static constexpr unsigned int toerrs = 300; //run until toerrs bit errors occur
static constexpr double amplitude = 1.0;

///run simulation with perfect channel knowledge
void runperfectchannelsim(const function<CodedConstellation*()> codedconstf, const vector<double>& snrdbs, const string name) {

  cout << "Running simulation " << name << " " << flush;
  auto starttime = chrono::system_clock::now();

  vector<future<double>> pearray(snrdbs.size()); //store the mse
  for( int snrind = 0; snrind < snrdbs.size(); snrind++ ) {
    //run each snr in a separate thread
    pearray[snrind] = async( launch::async, [=] {
  	CodedConstellation* codec = codedconstf();
  	const unsigned int K = codec->K; //number of information bits
  	const unsigned int absD = codec->numSymbols(); //number of data symbols
  	const unsigned int M = codec->sizeOfConstellation(); //M-PSK assumed
  	const unsigned int L = absD; //only data symbols in this case
	
  	std::vector<int> D;
  	for(int i = 0; i < L; i++) D.push_back(i); //data symbol positions (just 0 to L-1) in this case
    
  	double snrdB = snrdbs[snrind]; 
  	double var = amplitude*amplitude*pow(10, -snrdB/10); //variance of real and imaginary parts of noise
  	default_random_engine gen((long)clock()); //random number generator, seed with clock 
  	uniform_int_distribution<int> unifM(0, M-1); //for generating M-PSK symbols
  	uniform_real_distribution<double> unifphase(-pi,pi); //generator for true timeoffset 
  	normal_distribution<double> gn(0.0,sqrt(var));

  	vector<unsigned int> info(K); //vector of information bits
  	vector<complexd> s(L); //vector of symbols
  	vector<complexd> y(L); //vector of received symbols

  	//construct the channel inverter and decoder
  	InvertAndDecode invdec(codec, D);
    
  	int numerrs = 0; //count number of bit errors
  	int itrcount = 0; //count total number of iterations
  	while( numerrs < toerrs ) {
  	  const complex<double> a0 = polar<double>(amplitude, unifphase(gen)); //amplitude random phase
  	  for(int i = 0; i < K; i++) info[i] = (unsigned int) unifM(gen); //generate some random bits
  	  const vector<complexd>& cw = codec->encode(info); //encode the bits, result goes in cw
  	  for(int i = 0; i < L; i++) s[D[i]] = cw[i]; //convert to M-psk symbols
  	  for(int i = 0; i < L; i++) y[i] = a0*s[i] + complexd(gn(gen), gn(gen)); //compute transmitted signal
  	  const std::vector<unsigned int>& rbits = invdec.decode(y, a0, var); //decode recieved signal
  	  for(int i = 0; i < K; i++) numerrs += (rbits[i] == info[i]) ? 0 : 1;
  	  itrcount++;
  	}
	
  	delete codec;
  	cout << "." << flush;
  	return ((double)numerrs)/K/itrcount; 
      } );
  }

  //now write output to file
  for(int i = 0; i < pearray.size(); i++) pearray[i].wait();
  ofstream file(string("data/") + name);      
  for( int snrind = 0; snrind < snrdbs.size(); snrind++ ) 
    file << snrdbs[snrind] << "\t" << pearray[snrind].get() << endl;
  file.close();

  auto endtime = chrono::system_clock::now();
  auto elapsed_seconds = std::chrono::duration_cast<seconds>(endtime-starttime).count();
  cout << " finished in  " << elapsed_seconds << " seconds." << endl;

}

//run simulation with given initial channel estimator
void runsim(function<CodedConstellation*()> codedconstf, function<CoherentMackenthun(vector<int>&,vector<int>&,vector<complexd>&)> phestf, 
	    const vector<double>& snrdbs, const int numpilots, const string name) {
    
  cout << "Running simulation " << name << " " << flush;
  auto starttime = chrono::system_clock::now();
  
  vector<future<double>> pearray(snrdbs.size()); //store the mse
  for( int snrind = 0; snrind < snrdbs.size(); snrind++ ) {
    //run each snr in a separate thread
    pearray[snrind] = async( launch::async, [=] {
	
	CodedConstellation* codec = codedconstf();
  	const unsigned int K = codec->K; //number of information bits
  	const unsigned int absD = codec->numSymbols(); //number of data symbols
  	const unsigned int M = codec->sizeOfConstellation(); //M-PSK assumed
  	const unsigned int L = absD + numpilots; //total number of symbols
	
	std::vector<int> P;
	for(int i = 0; i < numpilots; i++) P.push_back(i); //pilots at the front
	std::vector<int> D;
	for(int i = numpilots; i < L; i++) D.push_back(i); //data at the back
	
	default_random_engine gen((long)clock());
	uniform_int_distribution<int> unifM(0, M-1); //for generating M-PSK symbols
	vector<complexd> pilots(numpilots); //vector of pilot symbols
	for(int i = 0; i < numpilots; i++) pilots[i] = polar<double>(1.0, (2*pi*unifM(gen)) / M);
	
	//construct phase estimator
	CoherentMackenthun phest = phestf(D,P,pilots);
	
	double snrdB = snrdbs[snrind]; 
	double var = amplitude*amplitude*pow(10, -snrdB/10); //variance of real and imaginary parts of noise
	uniform_real_distribution<double> unifphase(-pi,pi); //generator for true timeoffset 
	normal_distribution<double> gn(0.0,sqrt(var));

	//memory and pointers for LDPC bits
	vector<unsigned int> info(K); //vector of information bits
	vector<complexd> s(L); //vector of symbols
	vector<complexd> y(L); //vector of received symbols	
	for(int i = 0; i < numpilots; i++) s[P[i]] = pilots[i]; //fill pilots

	//construct the turbo synchroniser
	TurboSyncroniser turbosync(codec, D, P, pilots);
	
	int countcodeworderrors = 0;
	int numerrs = 0; //count number of bit errors
	int itrcount = 0; //count total number of iterations
	while( numerrs < toerrs ) {
	  const complex<double> a0 = polar<double>(amplitude, unifphase(gen)); //amplitude random phase
	  for(int i = 0; i < K; i++) info[i] = (unsigned int) unifM(gen); //generate some random bits
	  const vector<complexd>& cw = codec->encode(info); //encode the bits, result goes in cw
	  for(int i = 0; i < absD; i++) s[D[i]] = cw[i]; //convert to M-psk symbols
	  for(int i = 0; i < L; i++) y[i] = a0*s[i] + complexd(gn(gen), gn(gen)); //compute transmitted signal
	  phest.estimate(y); //run channel estimator
	  complexd ahat = phest.complexGainEstimate();
	  double varhat = phest.noiseVarianceEstimate();
	  const std::vector<unsigned int>& rbits = turbosync.decode(y, ahat, varhat); //decode recieved signal
	  for(int i = 0; i < K; i++) numerrs += (rbits[i] == info[i]) ? 0 : 1;
	  itrcount++;
	}

	delete codec;
	cout << "." << flush;
	return ((double)numerrs)/K/itrcount;
      } );
  }

  //now write output to file
  for(int i = 0; i < pearray.size(); i++) pearray[i].wait();
  ofstream file(string("data/") + name);      
  for( int snrind = 0; snrind < snrdbs.size(); snrind++ ) 
    file << snrdbs[snrind] << "\t" << pearray[snrind].get() << endl;
  file.close();

  auto endtime = chrono::system_clock::now();
  auto elapsed_seconds = std::chrono::duration_cast<seconds>(endtime-starttime).count();
  cout << " finished in  " << elapsed_seconds << " seconds." << endl;

}

int main(int argc, char** argv) {
    
  int M = 2; //QPSK at the moment

  vector<double> snrdbs;//snrs in db we will run
  for(double db = 0; db <= 3.76; db+=0.25) snrdbs.push_back(db);
  
  auto perfectchan = [&] () { return new CodedBPSK("C/RA1N128.dec"); };

  //run simulation with perfect channel knowledge
  runperfectchannelsim(perfectchan, snrdbs, "perfectchannel");
  
  //run turbo synchroniser initialised with least squares estimator
  auto cmack = [&] (vector<int>& D,vector<int>& P,vector<complexd>& p) { return CoherentMackenthun(D,P,p,M); };
  runsim(perfectchan, cmack, snrdbs, 5, "cmack5");
  //runsim("C/RA1N128.dec", cmack, snrdbs, 5, "cmack5");
  //runsim("C/RA1N128.dec", cmack, snrdbs, 20, "cmack20");
  //runsim("C/RA1N128.dec", cmack, snrdbs, 5, "cmack40");

  //run turbo synchroniser initial with pilots only
  auto pilotonly = [&] (vector<int>& D,vector<int>& P,vector<complexd>& p) { return CoherentMackenthun(vector<int>(),P,p,M); };
  runsim(perfectchan, pilotonly, snrdbs, 5, "pilotonly5");
  //runsim("C/RA1N128.dec", pilotonly, snrdbs, 5, "pilotonly5");
  //runsim("C/RA1N128.dec", pilotonly, snrdbs, 20, "pilotonly20");  
  //runsim("C/RA1N128.dec", pilotonly, snrdbs, 40, "pilotonly40");
}
