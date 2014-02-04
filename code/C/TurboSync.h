/* 
 * File:   TurboSync.h
 * Author: Robby McKilliam
 *
 * Created on 30 January 2014, 10:51 AM
 */

#ifndef TURBOSYNC_H
#define	TURBOSYNC_H

#include "CodedConstellation.h"
#include "CoherentMackenthun.h"
#include "LDPCDec.h"

using namespace std;

/** Base class for coded BPSK transmission */
class CodedBandpassReciever {
    
public:
    
    /// maximum number of decoder iterations.
    const unsigned int maxiterations;
    
    CodedBandpassReciever() = delete;
    
    /**
     * Decode and return recieved bits from recieved signal y given channel estimate chat and
     * noise variance estimate varhat.  varhat is the the variance of the complex noise, i.e., the
     * var(re) + var(im).
     */
    virtual const vector<unsigned int>& decode(const vector<complexd>& y, const complexd chat, const double varhat) = 0;
    
protected:
    
     CodedConstellation* codec;
    
     CodedBandpassReciever(CodedConstellation* decoder, const unsigned int maxitr=100) : codec(decoder), maxiterations(maxitr) {}
    
};

/** Standard channel inversion and decoding. Uses and LDPC code and assumes BPSK */
class InvertAndDecode : public CodedBandpassReciever {
public:

    InvertAndDecode(CodedConstellation* decoder, const std::vector<int>& Din, const unsigned int maxitr = 100) :
    CodedBandpassReciever(decoder, maxitr),
    D(Din),
    r(Din.size()) {
        if (D.size() != decoder->numSymbols()) throw "The number of data symbols must be the same as the number of symbols in a codeword";
    }

    virtual const std::vector<unsigned int>& decode(const std::vector<complexd>& y, const complexd chat, const double varhat) {
        for (int i = 0; i < D.size(); i++) r[i] = y[D[i]] / chat;
        return codec->decode(r, varhat/2, maxiterations); //run decoder, divide by 2 since varhat is the variance of the complex noise.
    }

protected:

    const std::vector<int> D; //data positions
    std::vector<complexd> r; //holds received signal after channel is inverted

};

/** Implements a BPSK turbo synchroniser using an LDPC code for phase and amplitude.  Assumes BPSK.  */
class TurboSyncroniser : public CodedBandpassReciever {
public:
        
    /** 
     * Construct a turbo synchroniser  using a LDPC decoder
     * and data symbols in positions described by the set D
     */
    TurboSyncroniser(CodedConstellation* decoder, const std::vector<int>& Din, const std::vector<int>& Pin, const std::vector<complexd> pilotsin, const unsigned int maxitr=100, const unsigned int turboitr=4) : 
    CodedBandpassReciever(decoder,maxitr), D(Din), P(Pin), pilots(pilotsin), L(Din.size()+Pin.size()), turboiterations(turboitr), ydata(Din.size()), fakeret(decoder->K) { }

    virtual const std::vector<unsigned int>& decode(const std::vector<complexd>& y, const complexd chat, const double varhat) {
        A = computeA(y);
        Ypilots = computeYpilots(y);
        //cout << A << ", " << Ypilots << ", " << Ypilots/((double)P.size()) << ", " << abs(Ypilots) << endl;
        return decoderrec(y, chat, varhat/2, turboiterations); //divide by 2 since varhat is the variance of the complex noise.
    }
    
     ///total number of transmitted symbols
    const int L;
    
    //number of turbo iteration to perform
    const unsigned int turboiterations;

protected:

    const std::vector<int> D; //data positions
    const std::vector<complexd> pilots; //pilots
    const std::vector<int> P; //data positions
    complexd Ypilots; //stores the sum of y corresponding with pilots
    double A; //store norm of received signal
    std::vector<complexd> ydata; //memory for data part of received signal
    std::vector<unsigned int> fakeret;

    /** Recursively run the turbo decoder */
    const std::vector<unsigned int>& decoderrec(const std::vector<complexd>& y, const complexd chat, const double varhat, const int itrcount){
        for(int i = 0; i < D.size(); i++) ydata[i] = y[D[i]]/chat; //data symbols with channel removed (by estimate))
        if(itrcount==0) {
            const std::vector<unsigned int>& r = codec->decode(ydata, varhat, maxiterations);
            for(int i = 0; i < r.size(); i++) fakeret[i] = r[i]; //copy memory, I have no idea why this is necessary, since the memory in codec should be just as good.  However, without this I was getting segfaults for reasons I can't explain!
            return fakeret;
        }
        //otherwise update channel and iterate
        complexd Y = Ypilots;
        const vector<complexd>& s = codec->expected(ydata,varhat,maxiterations); //expected symbols (soft symbols)
        for(int i = 0; i < D.size(); i++)  Y += y[D[i]]*conj(s[i]);
        complexd chatnew = Y/((double)L);
        double varhatnew = (A - std::norm(Y)/L) / L;
        decoderrec(y, chatnew, varhatnew, itrcount-1);
    }
    
    complexd computeYpilots(const std::vector<complexd>& y) const {
        complexd ret(0,0);
        for(int i = 0; i < P.size(); i++) ret += y[P[i]]*conj(pilots[i]);
        return ret;
    }
    
    double computeA(const std::vector<complexd>& y) const {
        double ret = 0.0; 
        for(int i = 0; i < P.size(); i++) ret += std::norm(y[P[i]]);
        for(int i = 0; i < D.size(); i++) ret += std::norm(y[D[i]]);
        return ret;
    }
    
};


#endif	/* TURBOSYNC_H */

