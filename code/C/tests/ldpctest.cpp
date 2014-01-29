/* 
 * File:   ldpctest.cpp
 * Author: Robby McKilliam
 *
 * Created on 29/01/2014, 1:58:37 PM
 */

#include <stdlib.h>
#include <iostream>
#include <functional>
#include <random>
#include "LDPCDec.h"

#define	TOL 1e-7

using namespace std;

/*
 * Test Gottfried's LDPC code
 */

bool testconstruct() {
    CLDPCDec codec = CLDPCDec("RA1N128.dec");
    //cout << codec.getN() << ", " << codec.getK() << ", " << codec.getE() << ", " << codec.getM() << endl;
    return codec.dec_loaded();
}

bool testencode() {
    CLDPCDec codec = CLDPCDec("RA1N128.dec");
    unsigned int *codeword = (unsigned int*) malloc(codec.getN() * sizeof (unsigned int));
    unsigned int *info = codeword;
    unsigned int *parity = codeword+codec.getK();
    for(int i = 0; i < codec.getK(); i++) info[i] = rand()%2;
    codec.encodeRA(info,parity);
    bool pass = true;
    for(int i = 0; i < codec.getK(); i++) pass &= codeword[i] == info[i];
    return pass;
}

bool testencodedecode() {
    CLDPCDec codec = CLDPCDec("RA1N128.dec");
    unsigned int *codeword = (unsigned int*) malloc(codec.getN() * sizeof (unsigned int));
    unsigned int *info = codeword;
    unsigned int *parity = codeword+codec.getK();
    for(int i = 0; i < codec.getK(); i++) info[i] = rand()%2;
    codec.encodeRA(info,parity);
    
    double *Lch = (double*) malloc(codec.getN() * sizeof (double));
    double *Lapp = (double*) malloc(codec.getN() * sizeof (double));
    
    //make the Lch super confident
    for(int i = 0; i < codec.getN(); i++) 
        if(codeword[i]==0) Lch[i] = 1000.0; else Lch[i] = -1000.0;
    
    codec.decode(Lch,Lapp,10); //run at most 10 iterations
    
    //check that Lapp makes sense
    bool pass = true;
    for(int i = 0; i < codec.getN(); i++) 
        pass &= ((codeword[i]==0)&&(Lapp[i] > 0)) || ((codeword[i]==1)&&(Lapp[i] < 0));
    
    return pass;
}

bool testencodedecodeAWGN() {
    double snrdB = 10.0;
    double amplitude = 1.0; //power of BPSK transmission
    double var = amplitude*amplitude*pow(10, -snrdB/10);
    
    //get code
    CLDPCDec codec = CLDPCDec("RA1N128.dec");
    
    //encode
    unsigned int *codeword = (unsigned int*) malloc(codec.getN() * sizeof (unsigned int));
    unsigned int *info = codeword;
    unsigned int *parity = codeword+codec.getK();
    for(int i = 0; i < codec.getK(); i++) info[i] = rand()%2;
    codec.encodeRA(info,parity);
    
    //add noise to get received signal
    default_random_engine generator;
    normal_distribution<double> distribution(0.0,sqrt(var));
    double *y = (double*) malloc(codec.getN() * sizeof (double));
    for(int i = 0; i < codec.getN(); i++) y[i] = 2.0*amplitude*(0.5-codeword[i]) + distribution(generator);
    
    //memory for llrs and fill channel llrs
    double *Lch = (double*) malloc(codec.getN() * sizeof (double));
    double *Lapp = (double*) malloc(codec.getN() * sizeof (double));
    for(int i = 0; i < codec.getN(); i++) Lch[i] = CLDPCDec::llrBPSK(y[i],amplitude,var);
    
    //run decoder for 50 iterations
    codec.decode(Lch,Lapp,50);
    
    bool pass = true;
    for(int i = 0; i < codec.getN(); i++) 
        pass &= ((codeword[i]==0)&&(Lapp[i] > 0)) || ((codeword[i]==1)&&(Lapp[i] < 0));
    
    return pass;
}

void runtest(string name, function<bool() > test) {
    cout << name << " ... ";
    if (!test()) cout << "FAIL" << endl;
    else cout << "pass" << endl;
}

int main(int argc, char** argv) {
    runtest("test construct LDPC", testconstruct);
    runtest("test encode RA", testencode);
    runtest("test encode and decode RA", testencodedecode);
    runtest("test encode and decode RA in AWGN", testencodedecodeAWGN);
    return (EXIT_SUCCESS);
}

