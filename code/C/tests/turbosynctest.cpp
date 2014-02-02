/* 
 * File:   turbosynctest.cpp
 * Author: Robby McKilliam
 *
 * Created on 30/01/2014, 4:46:30 PM
 */

#include <stdlib.h>
#include <iostream>
#include <functional>
#include <random>
#include "TurboSync.h"
#include "CoherentMackenthun.h"
#include "CodedConstellation.h"

#define	TOL 1e-7

using namespace std;

bool testInvertAndDecode() {
    
    CodedBPSK codec("RA1N128.dec");

    int M = 2;
    int absP = 50; //number of pilot symbols
    int absD = codec.N; //number of data symbols.  BPSK so one symbol per bit
    int L = absP + absD; //total number of symbols
    complexd a0 = polar<double>(1.0, 0.6); //the complex amplitude
    std::vector<int> P;
    for(int i = 0; i < absP; i++) P.push_back(i); //pilots at the front
    std::vector<int> D;
    for(int i = absP; i < L; i++) D.push_back(i); //data at the back

    default_random_engine generator;
    uniform_int_distribution<int> unifM(0, M-1); //for generating M-PSK symbols
    vector<complexd> pilots(absP); //vector of pilot symbols
    for(int i = 0; i < absP; i++) pilots[i] = polar<double>(1.0, 2 * pi * unifM(generator) / M);
    
    vector<unsigned int> info(codec.K);
    for(int i = 0; i < codec.K; i++) info[i] = (unsigned int) unifM(generator); //generate some random bits
    const vector<complexd>& codeword = codec.encode(info); //encode the bits, result goes in codeword
    
    vector<complexd> s(L); //vector of symbols
    for(int i = 0; i < absP; i++) s[P[i]] = pilots[i]; //pilots
    for(int i = 0; i < absD; i++) s[D[i]] = codeword[i]; //codeword
    
    double var0 = 0.001;
    normal_distribution<double> gn(0.0,sqrt(var0)); //Gaussian noise
    std::vector<complexd> y(L);
    for(int i = 0; i < L; i++) y[i] = a0*s[i] + complexd(gn(generator), gn(generator)); //compute transmitted signal
    
    InvertAndDecode invdec(&codec, D);
    const std::vector<unsigned int>& bits = invdec.decode(y,a0,var0);
    
    bool pass = true;
    for(int i = 0; i < codec.K; i++) pass &= (bits[i] == info[i]); //check decoded correctly
   
    return pass;
}

bool testTurboSyncWithPerfectInitialiser() {
    
    CodedBPSK codec("RA1N128.dec");

    int M = 2;
    int absP = 50; //number of pilot symbols
    int absD = codec.N; //number of data symbols.  BPSK so one symbol per bit
    int L = absP + absD; //total number of symbols
    complexd a0 = polar<double>(1.0, 0.6); //the complex amplitude
    std::vector<int> P;
    for(int i = 0; i < absP; i++) P.push_back(i); //pilots at the front
    std::vector<int> D;
    for(int i = absP; i < L; i++) D.push_back(i); //data at the back

    default_random_engine generator;
    uniform_int_distribution<int> unifM(0, M-1); //for generating M-PSK symbols
    vector<complexd> pilots(absP); //vector of pilot symbols
    for(int i = 0; i < absP; i++) pilots[i] = polar<double>(1.0, 2 * pi * unifM(generator) / M);
    
    vector<unsigned int> info(codec.K);
    for(int i = 0; i < codec.K; i++) info[i] = (unsigned int) unifM(generator); //generate some random bits
    const vector<complexd>& codeword = codec.encode(info); //encode the bits, result goes in codeword
    
    vector<complexd> s(L); //vector of symbols
    for(int i = 0; i < absP; i++) s[P[i]] = pilots[i]; //pilots
    for(int i = 0; i < absD; i++) s[D[i]] = codeword[i]; //codeword
    
    double var0 = 0.001;
    normal_distribution<double> gn(0.0,sqrt(var0)); //Gaussian noise
    std::vector<complexd> y(L);
    for(int i = 0; i < L; i++) y[i] = a0*s[i] + complexd(gn(generator), gn(generator)); //compute transmitted signal
    
    TurboSyncroniser turbosync(&codec, D, P, pilots);
    const std::vector<unsigned int>& bits = turbosync.decode(y,a0,var0);
    
    bool pass = true;
    for(int i = 0; i < codec.K; i++) pass &= (bits[i] == info[i]); //check decoded correctly
   
    return pass;
    
}

void runtest(string name, function<bool() > test) {
    cout << name << " ... ";
    if (!test()) cout << "FAIL" << endl;
    else cout << "pass" << endl;
}

int main(int argc, char** argv) {
    runtest("test invert and decode", testInvertAndDecode);
    runtest("test turbo synchroniser with perfect channel initialiser", testTurboSyncWithPerfectInitialiser);
    return (EXIT_SUCCESS);
}
