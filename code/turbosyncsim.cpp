/* 
 * File:   turbosyncsim.cpp
 * Author: Robby McKilliam
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "CoherentMackenthun.h"
#include "LDPCDec.h"
#include "TurboSync.h"
#include <random>
#include <functional>
#include <future>
#include <chrono>

using namespace std;

typedef std::chrono::duration<double, std::ratio<1>> seconds;

static constexpr double pi = 3.141592653589793238463;

static constexpr unsigned int M = 2; //BPSK (you can't actually change this!)
static constexpr unsigned int toerrs = 5000; //run until toerrs bit errors occur
static constexpr double absD = 128;
static constexpr double amplitude = 1.0;

///run simulation with perfect channel knowledge
void runperfectchannelsim(string decoderfilename, const vector<double>& snrdbs, string name) {

  cout << "Running simulation " << name << " " << flush;
  auto starttime = chrono::system_clock::now();

  vector<future<double>> pearray(snrdbs.size()); //store the mse
  for( int snrind = 0; snrind < snrdbs.size(); snrind++ ) {
    //run each snr in a separate thread
    pearray[snrind] = async( launch::async, [=] {
	CLDPCDec codec(decoderfilename.c_str());
	const unsigned int L = codec.getN();
	const unsigned int K = codec.getK();

	std::vector<int> D;
	for(int i = 0; i < L; i++) D.push_back(i); //data at the back
    
	double snrdB = snrdbs[snrind]; 
	double var = amplitude*amplitude*pow(10, -snrdB/10); //variance of real and imaginary parts of noise
	default_random_engine gen((long)clock()); //random number generator, seed with clock 
	uniform_int_distribution<int> unifM(0, M-1); //for generating M-PSK symbols
	uniform_real_distribution<double> unifphase(-pi,pi); //generator for true timeoffset 
	normal_distribution<double> gn(0.0,sqrt(var));

	//cout << var << endl;

	//memory and pointers for LDPC bits
	unsigned int *codeword = (unsigned int*) malloc(L * sizeof (unsigned int));
	unsigned int *info = codeword;
	unsigned int *parity = codeword+K;
	vector<complexd> s(L); //vector of symbols
	vector<complexd> y(L); //vector of received symbols

	//construct the channel inverter
	InvertAndDecode invdec(&codec, D);
    
	int countcodeworderrors = 0;
	int numerrs = 0; //count number of bit errors
	int itrcount = 0; //count total number of iterations
	while( numerrs < toerrs ) {
	  const complex<double> a0 = polar<double>(1.0, unifphase(gen)); //amplitude random phase
	  for(int i = 0; i < K; i++) info[i] = (unsigned int) unifM(gen); //generate some random bits
	  codec.encodeRA(info, parity); //encode the bits, result goes in codeword
	  for(int i = 0; i < L; i++) s[D[i]] = polar<double>(1.0, (2*pi*codeword[i]) / M); //convert to M-psk symbols
	  for(int i = 0; i < L; i++) y[i] = a0*s[i] + complexd(gn(gen), gn(gen)); //compute transmitted signal
	  const std::vector<unsigned int>& rbits = invdec.decode(y, a0, var); //decode recieved signal
	  for(int i = 0; i < K; i++) numerrs += (rbits[i] == info[i]) ? 0 : 1;
	  itrcount++;
	}
	//cout << snrdB << ", " << pearray[snrind] << endl;
	cout << "." << flush;
	return ((double)numerrs)/K/itrcount; } );
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
    
  vector<double> snrdbs;//snrs in db we will run
  for(double db = 0; db <= 4; db+=0.25) snrdbs.push_back(db);
  
  runperfectchannelsim("C/RA1N128.dec", snrdbs, "perfectchannel");
    
    
}
