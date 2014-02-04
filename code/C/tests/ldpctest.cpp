/* 
 * File:   ldpctest.cpp
 * Author: Robby McKilliam
 *
 * Created on 29/01/2014, 1:58:37 PM
 */

#include "LDPCDec.h"
#include "CodedConstellation.h"
#include <stdlib.h>
#include <iostream>
#include <functional>
#include <random>

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
    for(int i = 0; i < codec.getN(); i++) Lch[i] = CLDPCDec::BPSK2LLR(y[i],var);
    
    //run decoder for 50 iterations
    codec.decode(Lch,Lapp,50);
    
    bool pass = true;
    for(int i = 0; i < codec.getN(); i++) 
        pass &= ((codeword[i]==0)&&(Lapp[i] > 0)) || ((codeword[i]==1)&&(Lapp[i] < 0));
    
    return pass;
}

bool testCodedBPSKConstruct() {
    CodedBPSK cbpsk("RA1N128.dec");
    return true;
}

bool testCodedBPSKEncode() {
    CodedBPSK cbpsk("RA1N128.dec");
    vector<unsigned int> infobits(cbpsk.K);
    for(int i = 0; i < cbpsk.K; i++) infobits[i] = rand()%2;
    const vector<complexd>& cw = cbpsk.encode(infobits);
    bool pass = true;
    for(int i = 0; i < cbpsk.K; i++) {
        bool azero = (abs(cw[i]-complexd(1,0))<TOL) && (infobits[i] == 0);
        bool aone = (abs(cw[i]-complexd(-1,0))<TOL) && (infobits[i] == 1);
        pass &= azero || aone;
    }
    return pass;
}

bool testCodedBPSKAWGN() {
    double snrdB = 10.0;
    double amplitude = 1.0; //power of BPSK transmission
    double var = amplitude*amplitude*pow(10, -snrdB/10);
    default_random_engine gen;
    normal_distribution<double> gn(0.0,sqrt(var));
    
    CodedBPSK cbpsk("RA1N128.dec");
    vector<unsigned int> infobits(cbpsk.K);
    for(int i = 0; i < cbpsk.K; i++) infobits[i] = rand()%2;
    const vector<complexd>& cw = cbpsk.encode(infobits);
    
    //generate received signal with noise
    vector<complexd> y(cbpsk.numSymbols());
    for(int i = 0; i < cbpsk.numSymbols(); i++) y[i] = cw[i] + complexd(gn(gen),gn(gen));
    
    //decode signal
    const vector<unsigned int>& decodedbits = cbpsk.decode(y,var);
    
    bool pass = true;
    for(int i = 0; i < cbpsk.K; i++)
        pass &= decodedbits[i] == infobits[i];
    
    return pass;
}

/** Class for extracting protected methods for testing */
class TestCodedBPSK : public CodedBPSK {
public:
    TestCodedBPSK(const char* ldpcspec) : CodedBPSK(ldpcspec) {}
    const vector<complexd>& testLLRs2constellation(double* llrs) {
        return LLRs2constellation(llrs);
    }
};

bool testCodedBPSKLLR2Expected() {
    const int N = 256;
    vector<double> llr(N);
    llr[0] = -10000.0;    
    llr[1] = 10000.0;
    llr[2] = 0.0;
    TestCodedBPSK cbpsk("RA1N128.dec");
    auto s = cbpsk.testLLRs2constellation(&llr[0]);
    //cout << s[0] << ", " << s[1] << ", " << s[2] << endl;
    bool pass = std::abs(s[0]-complexd(-1.0,0)) < 0.01;
    pass &= std::abs(s[1]-complexd(1.0,0)) < 0.01;
    pass &= std::abs(s[2]-complexd(0.0,0)) < 0.01;
    return pass;
}

/** Class for extracting protected methods for testing */
class TestCodedQPSK : public CodedQPSK {
public:
    TestCodedQPSK(const char* ldpcspec) : CodedQPSK(ldpcspec) {}
    const vector<complexd>& testLLRs2constellation(double* llrs) {
        return LLRs2constellation(llrs);
    }
    const vector<complexd>& testcodewordbits2constellation(unsigned int* cw) {
        return codewordbits2constellation(cw);
    }
};

bool testCodedQPSKcodeworbits2constellation() {
    const int N = 256;
    vector<unsigned int> b(N);
    b[0]=0; b[1]=0; b[2]=0; b[3]=1; b[4]=1; b[5]=0; b[6]=1; b[7]=1;
    TestCodedQPSK cbpsk("RA1N128.dec");
    auto s = cbpsk.testcodewordbits2constellation(&b[0]);
    //cout << s[0] << ", " << s[1] << ", " << s[2] << endl;
    bool pass = std::abs(s[0]-complexd(1.0,1.0)/sqrt(2)) < TOL;
    pass &= std::abs(s[1]-complexd(1.0,-1.0)/sqrt(2)) < TOL;
    pass &= std::abs(s[2]-complexd(-1.0,1.0)/sqrt(2)) < TOL;
    pass &= std::abs(s[3]-complexd(-1.0,-1.0)/sqrt(2)) < TOL;
    return pass;
}

bool testCodedQPSKLLRs2constellation() {
    const int N = 256;
    const double g = 10000.0;
    vector<double> llr(N);
    llr[0]=g; llr[1]=g; llr[2]=g; llr[3]=-g; llr[4]=-g; llr[5]=g; llr[6]=-g; llr[7]=-g;
    TestCodedQPSK cbpsk("RA1N128.dec");
    auto s = cbpsk.testLLRs2constellation(&llr[0]);
    //cout << s[0] << ", " << s[1] << ", " << s[2] << endl;
    bool pass = std::abs(s[0]-complexd(1.0,1.0)/sqrt(2)) < 0.01;
    pass &= std::abs(s[1]-complexd(1.0,-1.0)/sqrt(2)) < 0.01;
    pass &= std::abs(s[2]-complexd(-1.0,1.0)/sqrt(2)) < 0.01;
    pass &= std::abs(s[3]-complexd(-1.0,-1.0)/sqrt(2)) < 0.01;
    pass &= std::abs(s[4]-complexd(0,0)/sqrt(2)) < 0.01;
    pass &= std::abs(s[5]-complexd(0,0)/sqrt(2)) < 0.01;
    return pass;
}

bool testCodedQPSKAWGN() {
    double snrdB = 10.0;
    double amplitude = 1.0; //power of BPSK transmission
    double var = amplitude*amplitude*pow(10, -snrdB/10);
    default_random_engine gen;
    normal_distribution<double> gn(0.0,sqrt(var));
    
    CodedQPSK cqpsk("RA1N128.dec");
    vector<unsigned int> infobits(cqpsk.K);
    for(int i = 0; i < cqpsk.K; i++) infobits[i] = rand()%2;
    const vector<complexd>& cw = cqpsk.encode(infobits);
    
    //generate received signal with noise
    vector<complexd> y(cqpsk.numSymbols());
    for(int i = 0; i < cqpsk.numSymbols(); i++) y[i] = cw[i] + complexd(gn(gen),gn(gen));
    
    //decode signal
    const vector<unsigned int>& decodedbits = cqpsk.decode(y,var);
    
    bool pass = true;
    for(int i = 0; i < cqpsk.K; i++)
        pass &= decodedbits[i] == infobits[i];
    
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
    runtest("test construct coded BPSK", testCodedBPSKConstruct);
    runtest("test coded BPSK encode", testCodedBPSKEncode);
    runtest("test coded BPSK in AWGN", testCodedBPSKAWGN);
    runtest("test mapping LLR to BPSK", testCodedBPSKLLR2Expected);
    runtest("test mapping bits to symbols for QPSK", testCodedQPSKcodeworbits2constellation);
    runtest("test mapping LLR to QPSK", testCodedQPSKLLRs2constellation);
    runtest("test coded QPSK in AWGN", testCodedQPSKAWGN);
    return (EXIT_SUCCESS);
}

