/* 
 * File:   newsimpletest.cpp
 * Author: Robby McKilliam
 *
 * Created on 29/01/2014, 1:38:53 PM
 */

#include <stdlib.h>
#include <iostream>
#include <functional>
#include <random>
#include "CoherentMackenthun.h"

#define	TOL 1e-7

using namespace std;

/*
 * Test Coherent Mackenthun code.
 */

bool testEstimate() {
    int M = 4; //QPSK
    std::vector<int> P;
    P.push_back(0);
    std::vector<int> D;
    D.push_back(1);
    D.push_back(2);
    D.push_back(3);
    D.push_back(4);
    std::vector<complexd> pilots;
    pilots.push_back(complexd(1, 0));
    std::vector<complexd> y;
    y.push_back(complexd(-0.15662206555348682, 0.998001546268412));
    y.push_back(complexd(-0.15571514521240135, 0.993262721949701));
    y.push_back(complexd(0.9856798266529281, 0.12784677325244354));
    y.push_back(complexd(0.14524463440340715, -0.9886039113420713));
    y.push_back(complexd(0.14550394304366623, -0.9934311310398819));
    CoherentMackenthun cmack(D, P, pilots, M);
    cmack.estimate(y);
    complexd ahat = cmack.complexGainEstimate();
    complexd e1(-0.14618651229308086, 0.9917958274505988);
    if (abs(e1 - ahat) < TOL) return true;
    else return false;
}

bool testEstimateNoiseVar() {
    int M = 2;
    int absP = 50;
    int absD = 100;
    int L = absP + absD;
    complexd a0 = polar<double>(1.0, 0.6); //the complex amplitude
    std::vector<int> P;
    for(int i = 0; i < absP; i++) P.push_back(i); //pilots at the front
    std::vector<int> D;
    for(int i = absP; i < L; i++) D.push_back(i); //data at the back
    
    default_random_engine generator;
    uniform_int_distribution<int> unifM(1, M); //for generating M-PSK symbols
    vector<complexd> pilots(absP); //vector of pilot symbols
    for(int i = 0; i < absP; i++) pilots[i] = polar<double>(1.0, 2 * pi * unifM(generator) / M);
    
    vector<complexd> s(L); //vector of symbols
    for(int i = 0; i < absP; i++) s[P[i]] = pilots[i];
    for(int i = 0; i < absD; i++) s[D[i]] = polar<double>(1.0, 2 * pi * unifM(generator) / M); //random M-psk symbols
    
    double var0 = 0.001;
    normal_distribution<double> gn(0.0,sqrt(var0)); //Gaussian noise
    std::vector<complexd> y(L);
    for(int i = 0; i < L; i++) y[i] = a0*s[i] + complexd(gn(generator), gn(generator)); //compute transmitted signal
    
    CoherentMackenthun cmack(D, P, pilots, M);
    cmack.estimate(y);
    bool pass = std::norm(a0 - cmack.complexGainEstimate()) < 0.01;
    pass &= abs(2*var0 - cmack.noiseVarianceEstimate()) < 0.01;
    return pass;
}

bool testWithOnlyPilots() {
    int M = 2;
    int absP = 50;
    int absD = 100;
    int L = absP + absD;
    complexd a0 = polar<double>(1.0, 0.6); //the complex amplitude
    std::vector<int> P;
    for(int i = 0; i < absP; i++) P.push_back(i); //pilots at the front
    std::vector<int> D;
    for(int i = absP; i < L; i++) D.push_back(i); //data at the back
    
    default_random_engine generator;
    uniform_int_distribution<int> unifM(1, M); //for generating M-PSK symbols
    vector<complexd> pilots(absP); //vector of pilot symbols
    for(int i = 0; i < absP; i++) pilots[i] = polar<double>(1.0, 2 * pi * unifM(generator) / M);
    
    vector<complexd> s(L); //vector of symbols
    for(int i = 0; i < absP; i++) s[P[i]] = pilots[i];
    for(int i = 0; i < absD; i++) s[D[i]] = polar<double>(1.0, 2 * pi * unifM(generator) / M); //random M-psk symbols
    
    double var0 = 0.001;
    normal_distribution<double> gn(0.0,sqrt(var0)); //Gaussian noise
    std::vector<complexd> y(L);
    for(int i = 0; i < L; i++) y[i] = a0*s[i] + complexd(gn(generator), gn(generator)); //compute transmitted signal
    
    CoherentMackenthun cmack(vector<int>(), P, pilots, M);
    cmack.estimate(y);
    bool pass = std::norm(a0 - cmack.complexGainEstimate()) < 0.01;
    pass &= abs(2*var0 - cmack.noiseVarianceEstimate()) < 0.01;
    return pass;
}

void runtest(string name, function<bool()> test) {
    cout << name << " ... ";
    if (!test()) cout << "FAIL" << endl;
    else cout << "pass" << endl;
}

int main(int argc, char** argv) {
    runtest("test coherent Mackenthun estimate", testEstimate);
    runtest("test noise variance estimator", testEstimateNoiseVar);
    runtest("test pilot only estimator", testWithOnlyPilots);
    return (EXIT_SUCCESS);
}

