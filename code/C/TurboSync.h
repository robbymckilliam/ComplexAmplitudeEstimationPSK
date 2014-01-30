/* 
 * File:   TurboSync.h
 * Author: Robby McKilliam
 *
 * Created on 30 January 2014, 10:51 AM
 */

#ifndef TURBOSYNC_H
#define	TURBOSYNC_H

#include "LDPCDec.h"

/** Implements a BPSK turbo synchroniser for phase and amplitude */
template <class PhaseEstimator> //P should be a FinitePulse
class TurboSyncroniser {
public:

    /** 
     * Construct a turbo synchroniser wth phase and amplitude  estimate phasest, LDPC decoder
     * decoder and data symbols in positions described by the set D
     */
    TurboSyncroniser(const PhaseEstimator& phaseest, const CLDPCDec& decoder, const std::vector<int>& Din) :
    phaseest(pest),
    D(Din),
    codec(decoder) {
        if(D.size() != codec.getN()) throw "The number of data symbols must be the same as the length of the code";
    }

    virtual const std::vector<unsigned int>& decode(const std::vector<complexd>& y) {
        pest.estimate(y); //estimate the channel
        complexd chat = pest.complexGainEstimate();
        double noisevarest = pest.noiseVarianceEstimate();
        
    }

protected:
    const PhaseEstimator pest;
    const std::vector<int> D; //data positions
    const CLDPCDec codec;
    double *Lch;
    double *Lapp;


};


#endif	/* TURBOSYNC_H */

