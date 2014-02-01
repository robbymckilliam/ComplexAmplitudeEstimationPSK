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

using namespace std;

class CodedConstellation {
public:
    
    ///the LDPC codec
    const CLDPCDec ldpcdecoder;
    ///number of information bits
    const unsigned int K;
    ///number bits in a codeword
    const unsigned int N;
    
    CodedConstellation(const string ldpcspec) : 
        ldpcdecoder(ldpcspec), 
        K(ldpcdecoder.getK()),
        N(ldpcdecoder.getN()),
        codewordbits(N),
        infobits(K),
        Lapp(N)
    {
    }
    
    /** 
     * Copy constructor forbidden since I don't want to write a copy constructor for Gottfried's
     * LDPC class.
     */
    CodedConstellation(const CodedConstellation& orig) = delete;
    
    ///Encodes bits to constellation points
    virtual const vector<complexd>& encode(const vector<unsigned int>& bits) {
        vector<unsigned int>& cw = infobits2codewordbits(bits);
        return codewordbits2constellation(cw);
    }
    
    ///Decodes a sequence of complex numbers to bits.  Requires input channel variance var.
    virtual const vector<unsigned int>& decode(const vector<complexd>& r, double var, unsigned int iters=100) {
        const vector<double>& lch = constellation2LLRs(r, var);
        ldpcdecoder.decode(&lch[0],&Lapp[0],iters);
        for(int i = 0; i < N; i++) codewordbits[i] = LLR2bit(Lapp[i]);
        return codewordbits2infobits(codewordbits);
    }
    
    /**
     * Return expected value of symbols (soft symbols) after a specified number of
     * iterations of the decoder.  Default is 100 iteration
     */
    virtual const vector<complexd>& expected(const vector<complexd>& r, unsigned int iters=100 ) {
        const vector<double>& lch = constellation2LLRs(r);
        ldpcdecoder.decode(&lch[0],&Lapp[0],iters);
        return LLRs2constellation(Lapp);
    }
    
    ///Maps LLRs to bits.  Default is positive llrs to 0 and negative llrs to 1
    virtual const unsigned int LLR2bit(double llr) { return (llr > 0) ? 0 : 1; }
    
protected:
    
    //memory for output (mutable!)
    vector<unsigned int> codewordbits;
    vector<unsigned int> infobits;
    vector<double> Lapp;
    
    ///Maps a sequence of constellation points to a sequence of log likelihood ratios
    virtual const vector<double>& constellation2LLRs(const vector<complexd>& r, double var) = 0;
    
    ////Maps LLRs to expected points on the complex plain (soft decisions)
    virtual const vector<complexd>& LLRs2constellation(const vector<double>& llrs) = 0;
    
    ///Maps codeword to it's specified constellation point
    virtual const vector<complexd>& codewordbits2constellation(const vector<unsigned int>& cw) = 0;
    
    /** 
     * Map codeword bits to info bits.  By default a systematic LDPC is uses, so we need only
     * take the first K bits from the codeword.
     */
    virtual const vector<unsigned int>& codewordbits2infobits(const vector<unsigned int>& cw) {
        for(int i = 0; i < K; i++) infobits[i] = cw[i];
        return infobits;
    }
    
    ///Map info bits to codeword, this is the LDPC encoder
    virtual const vector<unsigned int>& infobits2codewordbits(const vector<unsigned int>& ib) {
        ldpcdecoder.encodeRA(&ib[0],&codewordbits[K]);
        for(int i = 0; i < K; i++) codewordbits[i] = ib[i]; //copy info bits to front of codeword (a systematic code))
        return codewordbits;
    }

};

/** Class for low density parity check coded binary phase shift keying */
class CodedBPSK : public CodedConstellation {
    
public:

    CodedBPSK(const string ldpcspec) : 
        CodedConstellation(ldpcspec),
         Lch(N),
        codeword(N)
    {}
    
protected:
    
    //memory for output (mutable!)
    vector<complexd> codeword;
    vector<double> Lch; 
    
    virtual const vector<double>& constellation2LLRs(const vector<complexd>& r, double var) {
        for(int i = 0; i < N; i++) Lch[i] = 2*real(codeword[i])/var;
        return Lch;
    }
    
    virtual const vector<complexd>& LLRs2constellation(const vector<double>& llrs) {
        for(int i = 0; i < N; i++) codeword[i] = atan(llrs[i]/2);
        return codeword;
    }
    
    virtual const vector<complexd>& codewordbits2constellation(const vector<unsigned int>& cw) {
        for(int i = 0; i < N; i++) codeword[i] = (cw[i]==0) ? complexd(1,0) : complexd(-1,0); 
        return codeword;
    }
    
};

#endif	/* CODEDCONSTELLATION_H */

