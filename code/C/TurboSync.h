/* 
 * File:   TurboSync.h
 * Author: Robby McKilliam
 *
 * Created on 30 January 2014, 10:51 AM
 */

#ifndef TURBOSYNC_H
#define	TURBOSYNC_H

#include "CoherentMackenthun.h"
#include "LDPCDec.h"

using namespace std;

/** Base class for coded BPSK transmission */
class CodedTransmission {
    
public:
    
    /// maximum number of decoder iterations.
    const unsigned int maxiterations;
    
    CodedTransmission(CLDPCDec* decoder, const unsigned int maxitr=30) : codec(decoder), maxiterations(maxitr) {
        Lch = (double*) malloc(codec->getN() * sizeof (double));
        Lapp = (double*) malloc(codec->getN() * sizeof (double));
    }
    
    ~CodedTransmission(){
       free(Lch);
       free(Lapp);
   }
    
    /**
     * Decode and return recieved bits from recieved signal y given channel estimate chat and
     * noise variance estimate varhat.
     */
    virtual const vector<unsigned int>& decode(const vector<complexd>& y, const complexd chat, const double varhat) = 0;
    
protected:
     CLDPCDec* codec;
     double *Lch; //memory for Gottfrieds decoder
     double *Lapp;
    
};

/** Standard channel inversion and decoding. Uses and LDPC code and assumes BPSK */
class InvertAndDecode : CodedTransmission {
    
public:
    
    InvertAndDecode(CLDPCDec* decoder, const std::vector<int>& Din, const unsigned int maxitr=30) : 
    CodedTransmission(decoder, maxitr),
    D(Din),
    bits(decoder->getK())
    {
        if(D.size() != codec->getN()) throw "The number of data symbols must be the same as the length of the code";
    }
    
    virtual const std::vector<unsigned int>& decode(const std::vector<complexd>& y, const complexd chat, const double varhat) {
        double rhohat = std::abs(chat); //channel amplitude estimate
        
        //fill channel log likelihood ratios
        for(int i = 0; i < codec->getN(); i++) {       
            double r = std::real(y[D[i]] / chat) * rhohat;
            Lch[i] = CLDPCDec::llrBPSK(r, rhohat, varhat/2); //divide noise by 2 since taking variance of real part for BPSK
        }
        codec->decode(Lch,Lapp,maxiterations); //run decoder
        
        //convert output log likelihood ratios to bits
        for(int i = 0; i < codec->getK(); i++) bits[i] = (Lapp[i] > 0) ? 0 : 1; 
        
        return bits;
    }
    
protected:
    std::vector<unsigned int> bits; //decoded information bits.
    const std::vector<int> D; //data positions
    
};

///** Implements a BPSK turbo synchroniser using an LDPC code for phase and amplitude */
//class TurboSyncroniser : CodedTransmission {
//public:
//
//    ///maximum number of turbo iterations.
//    const unsigned int maxiterations;
//    
//    /** 
//     * Construct a turbo synchroniser  using a LDPC decoder
//     * and data symbols in positions described by the set D
//     */
//    TurboSyncroniser(const CLDPCDec& decoder, const std::vector<int>& Pin, const std::vector<int>& Din, const std::vector<complexd> pilots, const unsigned int maxitr=30) :
//    D(Din),
//    codec(decoder),
//    maxiterations(maxitr),
//    bits(decoder.getK())
//    {
//        if(D.size() != codec.getN()) throw "The number of data symbols must be the same as the length of the code";
//        Lch = (double*) malloc(codec.getN() * sizeof (double));
//        Lapp = (double*) malloc(codec.getN() * sizeof (double));
//    }
//    
//   ~TurboSyncroniser(){
//       free(Lch);
//       free(Lapp);
//   }
//
//protected:
//    const std::vector<unsigned int> bits; //decoded information bits.
//    const std::vector<complexd> pilots; //pilots
//    const std::vector<int> P; //data positions
//    const std::vector<int> D; //data positions
//    const CLDPCDec codec;
//    double *Lch;
//    double *Lapp;
//
//
//};


#endif	/* TURBOSYNC_H */

