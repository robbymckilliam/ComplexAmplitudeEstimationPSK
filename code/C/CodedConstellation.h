/* 
 * File:   CodedConstellation.h
 * 
 * Contains classes for mapping binary LDPC codes to constellations in the complex plain
 * 
 * Author: Robby McKilliam
 *
 * Created on 31 January 2014, 9:35 PM
 */

#ifndef CODEDCONSTELLATION_H
#define	CODEDCONSTELLATION_H

#include "Util.h"
#include "LDPCDec.h"
#include <complex>
#include <vector>
#include <stdlib.h>
#include <iostream>


using namespace std;

class CodedConstellation {

protected:
    ///the LDPC codec
    CLDPCDec ldpcdecoder;    

public:
        
    ///number of information bits
    const unsigned int K;
    ///number bits in a codeword
    const unsigned int N;
    
    CodedConstellation(const string ldpcspec) :
        ldpcdecoder((char*)ldpcspec.c_str()),
        K(ldpcdecoder.getK()),
        N(ldpcdecoder.getN()),
        infobitsout(K)
    {
        Lch = (double*) malloc(N * sizeof (double));
        Lapp = (double*) malloc(N * sizeof (double));
        codewordbits = (unsigned int*) malloc(N * sizeof (unsigned int));
        infobits = codewordbits; 
        paritybits = codewordbits + K;
        //std::cout << N << ", " << K << ", " << infobitsout.size() << ", " << ldpcdecoder.dec_loaded() << std::endl;
    }
    
    ~CodedConstellation(){
       free(Lch);
       free(Lapp);
   }
    
    /** 
     * Copy constructor forbidden since I don't want to write a copy constructor for Gottfried's
     * LDPC class.
     */
    CodedConstellation(const CodedConstellation& orig) = delete;
    
    ///Encodes bits to constellation points
    virtual const vector<complexd>& encode(const vector<unsigned int>& bits) {
        for(int i = 0; i < K; i++) infobits[i] = bits[i]; //copy input bits to infobits for encode
        ldpcdecoder.encodeRA(infobits,paritybits); //encode
        return codewordbits2constellation(codewordbits);
    }
    
    ///Decodes a sequence of complex numbers to bits.  Requires input channel variance var.
    virtual const vector<unsigned int>& decode(const vector<complexd>& r, double var, unsigned int iters=100) {
        constellation2LLRs(r, var);
        ldpcdecoder.decode(Lch,Lapp,iters);
        for(int i = 0; i < K; i++) infobitsout[i] = LLR2bit(Lapp[i]);
        return infobitsout;
    }
    
    /**
     * Return expected value of symbols (soft symbols) after a specified number of
     * iterations of the decoder.  Default is 100 iteration
     */
    virtual const vector<complexd>& expected(const vector<complexd>& r, double var, unsigned int iters=100 ) {
        constellation2LLRs(r, var);
        ldpcdecoder.decode(Lch,Lapp,iters);
        return LLRs2constellation(Lapp);
    }
    
    ///Maps LLRs to bits.  Default is positive llrs to 0 and negative llrs to 1
    virtual const unsigned int LLR2bit(double llr) { return (llr > 0) ? 0 : 1; }
    
    ///Return the number of symbols in the constellation
    virtual const unsigned int numSymbols() = 0;
    
    ///Return the number of points in the complex constellation
    virtual const unsigned int sizeOfConstellation() = 0;
    
protected:
    
    //memory for bit output
    vector<unsigned int> infobitsout;
    
    //memory for Gottfried's decoder
    double *Lapp;
    double *Lch;
    unsigned int *codewordbits;
    unsigned int *infobits;
    unsigned int *paritybits;
    
    /** 
     * Maps a sequence of constellation points to a sequence of log likelihood ratios.
     * Result goes in Lch memory.  var is the variance of the real (or imaginary) part of the noise.
     */
    virtual void constellation2LLRs(const vector<complexd>& r, double var) = 0;
    
    ////Maps LLRs to expected points on the complex plain (soft decisions)
    virtual const vector<complexd>& LLRs2constellation(double* llrs) = 0;
    
    ///Maps codeword to its specified constellation point
    virtual const vector<complexd>& codewordbits2constellation(unsigned int* cw) = 0;

};

/** Class for low density parity check coded binary phase shift keying */
class CodedBPSK : public CodedConstellation {
    
public:
    
    CodedBPSK(const char* ldpcspec) : 
        CodedConstellation(ldpcspec),
        codeword(N)
    {}
        
   ///Number of symbols is the same as the length of the code for BPSK, one symbol per bit
    virtual const unsigned int numSymbols() { return N; }
    
   ///2 symbols in the BPSK constellation
   virtual const unsigned int sizeOfConstellation() { return 2; }
    
protected:
    
    //memory for output (mutable!)
    vector<complexd> codeword;
    
    virtual void constellation2LLRs(const vector<complexd>& r, double var) {
        for(int i = 0; i < N; i++) Lch[i] = 2*real(r[i])/var;
    }
    
    virtual const vector<complexd>& LLRs2constellation(double* llrs) {
        for(int i = 0; i < N; i++) codeword[i] = atan(llrs[i]/2)*2/pi;
        return codeword;
    }
    
    virtual const vector<complexd>& codewordbits2constellation(unsigned int* cw) {
        for(int i = 0; i < N; i++) codeword[i] = (cw[i]==0) ? complexd(1,0) : complexd(-1,0); 
        return codeword;
    }
    
};

/** Class for LDPC codes quaternary phase shift keying */
class CodedQPSK : public CodedConstellation {
public:

    CodedQPSK(const char* ldpcspec) :
    CodedConstellation(ldpcspec),
    codeword(N/2) {
        if (N%2 != 0) throw "Number of bits in codeword must be even for QPSK";
    }
    
    ///Number of symbols is half the length of the code for QPSK, one symbol per 2 bits
    virtual const unsigned int numSymbols() { return N/2; }
    
   ///4 symbols in the QPSK constellation
   virtual const unsigned int sizeOfConstellation() { return 4; }

protected:

    //memory for output (mutable!)
    vector<complexd> codeword;
    
    virtual void constellation2LLRs(const vector<complexd>& r, double var) {
        for(int i = 0; i < N/2; i++) {
            Lch[2*i] = 2*real(r[i])/var; //even bits encoded by real part
            Lch[2*i+1] = 2*imag(r[i])/var; //odd bits encoded by imaginary par
        }
    }
    
    virtual const vector<complexd>& codewordbits2constellation(unsigned int* cw) {
        for(int i = 0; i < N/2; i++) {
            double re = (cw[2*i]==0) ? sqrt(2)/2 : -sqrt(2)/2; //sqrt(2)/2 so that symbols have magnitude 1
            double im = (cw[2*i+1]==0) ? sqrt(2)/2 : -sqrt(2)/2; 
            codeword[i] = complexd(re,im); 
        }
        return codeword;
    }
    
    virtual const vector<complexd>& LLRs2constellation(double* llrs) {
        for(int i = 0; i < N/2; i++) {
            double re = sqrt(2)/pi * atan(llrs[2*i]/2);
            double im = sqrt(2)/pi * atan(llrs[2*i+1]/2);
            codeword[i] = complexd(re,im); 
        }
        return codeword;
    }

};

/** 
 * Class for LDPC coded QPSK with constellation rotated by pi/4 so that it's either on the real or
 * imaginary axis.
 */
class CodedQPSKRotated : public CodedQPSK {
    
public:
     
    const complexd rotate;
    
    CodedQPSKRotated(const char* ldpcspec) : CodedQPSK(ldpcspec) , rotate(std::polar<double>(1.0,pi/4)) {}
    
protected:
    
    virtual void constellation2LLRs(const vector<complexd>& r, double var) {
        for(int i = 0; i < N/2; i++) {
            complexd s = r[i] / rotate;
            Lch[2*i] = 2*real(s)/var; //even bits encoded by real part
            Lch[2*i+1] = 2*imag(s)/var; //odd bits encoded by imaginary par
        }
    }
    
    virtual const vector<complexd>& codewordbits2constellation(unsigned int* cw) {
        for(int i = 0; i < N/2; i++) {
            double re = (cw[2*i]==0) ? sqrt(2)/2 : -sqrt(2)/2; //sqrt(2)/2 so that symbols have magnitude 1
            double im = (cw[2*i+1]==0) ? sqrt(2)/2 : -sqrt(2)/2; 
            codeword[i] = complexd(re,im) * rotate; 
        }
        return codeword;
    }
    
    virtual const vector<complexd>& LLRs2constellation(double* llrs) {
        for(int i = 0; i < N/2; i++) {
            double re = sqrt(2)/pi * atan(llrs[2*i]/2);
            double im = sqrt(2)/pi * atan(llrs[2*i+1]/2);
            codeword[i] = complexd(re,im) * rotate; 
        }
        return codeword;
    }
    
};

#endif	/* CODEDCONSTELLATION_H */

