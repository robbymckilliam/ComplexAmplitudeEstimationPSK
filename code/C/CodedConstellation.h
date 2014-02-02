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
     * Result goes in Lch memory.
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

#endif	/* CODEDCONSTELLATION_H */

