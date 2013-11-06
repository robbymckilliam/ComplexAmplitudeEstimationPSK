/* 
 * File:   OerderMeyerMassey.h
 * Author: Robby McKilliam
 *
 * Created on 4 February 2013, 8:22 PM
 */

#ifndef OERDERMEYERMASSEY_H
#define	OERDERMEYERMASSEY_H

#include "TimeOffsetEstimator.h"
#include "FinitePulse.h"
#include "MainLoop.h"
#include <complex>
#include <vector>
#include <algorithm>

#ifndef VALTYPE
#define VALTYPE VALTYPE
#endif
#ifndef MSGTYPE
#define MSGTYPE VALTYPE
#endif


namespace TimeOffset {


/** 
 *Symbol timing estimate based on the filter and square estimator from:
 * 
 * Oerder and Meyr. "Digital Filter and Square Timing Recovery" IEEE TRANSACTIONS ON COMMUNICATIONS, VOL 36. NO.5, MAY 1988
 * 
 * @params S is the indices of all the transmitted symbols 
 * @params g is the transmission pulse shape, outside the interval [gtmin, gtmax] the pulse g is zero.
 * @params Ts is the sample period
 * @params T is the symbol period
 */
template <class Pulse>
class OerderMeyr : public TimeOffsetEstimator {
public:

    OerderMeyr(const vector<int>& S,
            const Pulse& g,
            const VALTYPE T,
            const VALTYPE Ts,
            const VALTYPE taumin,
            const VALTYPE taumax) :
    g(g), T(T), Ts(Ts), taumin(taumin), taumax(taumax) {
        if(taumax <= taumin) throw "Maximum time offset taumax must be larger than minimum time offset taumin";
        Smax = *std::max_element(S.begin(), S.end());
        Smin = *std::min_element(S.begin(), S.end());
    }

    /**
     * Runs the Oerder Meyr estimator. This returns an estimator of the symbol period modulo T.
     */
    virtual VALTYPE estimate(const vector<complex<VALTYPE> >& r) {
        setupr(r);
        int kmin = (int) ceil(4 * (g.tmin() + taumin + Smin * T) / T);
        int kmax = (int) floor(4 * (g.tmax() + taumax + Smax * T) / T);
        //sum squared matched filtered output complex exponential
        //complex<VALTYPE> X(0, 0);
        //for (int k = kmin; k <= kmax; k++)
        //    X += polar<VALTYPE>((VALTYPE)1.0, -(2 * pi * k) / 4) * norm(dotr(k * T / 4));
        VALTYPE Xreal = 0; VALTYPE Ximag = 0;
        for (int k = kmin; k <= kmax; k++) {
            VALTYPE d = norm(dotr(k * T / 4));
            int mod = k%4;
            if(mod==0) Xreal += d;
            else if(mod==1) Ximag -= d;
            else if(mod==2) Xreal -=d;
            else Ximag += d;
        }
        complex<VALTYPE> X(Xreal,Ximag);
        return -arg<VALTYPE>(X) / 2.0 / pi*T;
    }

    /** Inner product between the received signal and pulse g */
    virtual complex<VALTYPE> dotr(VALTYPE tau) const {
        int A = (int) ceil((g.tmin() + tau) / Ts);
        int B = (int) floor((g.tmax() + tau) / Ts);
        complex<VALTYPE> sum = complex<VALTYPE>(0, 0);
        for (int n = A; n <= B; n++) sum += r(n) * g.pulse(n * Ts - tau);
        return sum;
    }

    virtual string name() const {
        return "OerderMeyr";
    }
    
    virtual void setupr(const vector<complex<VALTYPE> >& r) {
        rmem = r.data();
        rmemsize = r.size();
    }

protected:
    int Smin;
    int Smax;
    const VALTYPE T;
    const VALTYPE Ts;
    const VALTYPE taumin;
    const VALTYPE taumax;
    const Pulse g;
    /** Pointer to current received signal */
    const complex<VALTYPE>* rmem;
    /** size of current received signal */
    unsigned int rmemsize;
    
    /** The received sequence, zeros returned for unknown values */
    inline complex<VALTYPE> r(int n) const {
        if (n > 0 && n <= rmemsize) return rmem[n - 1];
        else return complex<VALTYPE>(0, 0);
    }

};

/**
 * Runs the Oerder Meyer estimator followed by Massey's frame synchronisation algorithm to allow for complex
 * symbols.  The original paper is:
 * 
 * JAMES L. MASSEY, "Optimum Frame Synchronization", IEEE TRANSACTIONS ON COMMUNICATIONS, VOL. COM-20, NO. 2, APRIL 1972.
 * 
 * This version does not store the rho_k sequence, a faster version is OerderMeyerAndMassey below
 */
template <class Pulse>
class OerderMeyerAndMasseyNoStoreRho : public OerderMeyr<Pulse> {
    
public:

    OerderMeyerAndMasseyNoStoreRho(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<VALTYPE> >& pilots,
            const Pulse& g,
            const VALTYPE T,
            const VALTYPE Ts,
            const VALTYPE taumin,
            const VALTYPE taumax) : OerderMeyr<Pulse>(concatenate(P, D), g, T, Ts, taumin, taumax), P(P), D(D), pilots(pilots) {
    }

    virtual VALTYPE estimate(const vector<complex<VALTYPE> >& r) {
        setupgam(OerderMeyr<Pulse>::estimate(r)); //get estimate modulo the symbol period (call superclass)
        int imin = (int) floor((this->taumin - gam) / this->T);
        int imax = (int) ceil((this->taumax - gam) / this->T);
        int ibest = imin;
        VALTYPE maxL = -1; //the objective is always positive, so this is fine.
        for (int i = imin; i <= imax; i++) {
            VALTYPE thisL = L(i);
            if (thisL > maxL) {
                maxL = thisL;
                ibest = i;
            }
        }
        return ibest * this->T + gam;
    }

    /** The Massey-like objective function */
    VALTYPE L(const int d) const {
        VALTYPE dsum = 0.0;
        for (int i  = 0; i < D.size(); i++) dsum += std::abs<VALTYPE>(rho(D[i] + d));
        std::complex<VALTYPE> psum = std::complex<VALTYPE>(0, 0);
        for (int i = 0; i < P.size(); i++) psum += rho(P[i] + d) * std::conj<VALTYPE>(pilots[i]);
        return dsum + std::abs<VALTYPE>(psum);
    }

protected:
    const vector<int> P;
    const vector<int> D;
    const vector<complex<VALTYPE> > pilots;
    /** Variable contains the symbol time estimate.  MUTABLE */
    VALTYPE gam;
    
    virtual void setupgam(VALTYPE gam) {
        this->gam = gam;
    }
    
    inline virtual complex<VALTYPE> rho(const int k) const {
        return dotr((k + 1) * this->T + gam);
    }

};


/**
 * Runs the Oerder Meyer estimator followed by Massey's frame synchronisation algorithm to allow for complex
 * symbols.  The original paper is:
 * 
 * JAMES L. MASSEY, "Optimum Frame Synchronization", IEEE TRANSACTIONS ON COMMUNICATIONS, VOL. COM-20, NO. 2, APRIL 1972.
 * 
 * This version does not store the rho_k sequence, a faster version is OerderMeyerAndMassey below
 */
template <class Pulse>
class OerderMeyerAndMassey : public OerderMeyerAndMasseyNoStoreRho<Pulse> {
    
public:

    OerderMeyerAndMassey(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<VALTYPE> >& pilots,
            const Pulse& g,
            const VALTYPE T,
            const VALTYPE Ts,
            const VALTYPE taumin,
            const VALTYPE taumax) : OerderMeyerAndMasseyNoStoreRho<Pulse>(P,D,pilots,g,T,Ts,taumin,taumax) {
        int imin = (int) floor((taumin - 0.5) / T);
        int imax = (int) ceil((taumax + 0.5) / T);
        int kmin = imin + this->Smin;
        int kmax = imax + this->Smax;
        rhostore.resize(kmax - kmin + 1); //rhostore will never need to be bigger that this
    }

    virtual string name() const {
        return "OerderMeyrAndMassey";
    }

protected:
    vector<complex<VALTYPE> > rhostore;
    int rhooffset;
    
    /** Override setgam to also store all the values of rho we will need */
    virtual void setupgam(VALTYPE gam) {
        this->gam = gam;
        int imin = (int) floor((this->taumin - gam) / this->T);
        int imax = (int) ceil((this->taumax - gam) / this->T);
        int kmin = imin + this->Smin;
        int kmax = imax + this->Smax;
        rhooffset = kmin;
        for(int k = kmin; k <= kmax; k++) 
            rhostore[k - rhooffset] = OerderMeyerAndMasseyNoStoreRho<Pulse>::rho(k); //fill rhostore
    }
    
    /** Override rho to take from the array of stored values */
    inline virtual complex<VALTYPE> rho(const int k) const {
        return rhostore[k - rhooffset];
    }
            
};

#ifndef DEFAULTTABLESIZE
#define DEFAULTTABLESIZE 10000
#endif

/**
 * Tabulated Oerder Meyer.  Store pulse in a lookup table for fast access.  Table is memory
 * aligned for extra speed.
 */
template <class Pulse>
class TabulatedOerderMeyerAndMassey : public OerderMeyerAndMassey<Pulse> {
    
public:

    TabulatedOerderMeyerAndMassey(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<VALTYPE> >& pilots,
            const Pulse& g,
            const VALTYPE T,
            const VALTYPE Ts,
            const VALTYPE taumin,
            const VALTYPE taumax,
            const unsigned int mintabs = DEFAULTTABLESIZE) : 
    OerderMeyerAndMassey<Pulse>(P,D,pilots,g,T,Ts,taumin,taumax) {
        tabls = (int) ceil(mintabs*Ts/T);
        stepwidth = Ts/tabls;
        for (VALTYPE t = g.tmin(); t <= g.tmax(); t += stepwidth) pulsetable.push_back(g.pulse(t));
    }
    
    /** Inner product between the received signal and pulse g */
    inline virtual complex<VALTYPE> dotr(VALTYPE tau) const {
        //setup loop counters DOUBLES ARE DELIBERATE, DO NOT MODIFY
        int A = (int) ceil((this->g.tmin() + tau) / this->Ts);
        int B = (int) floor((this->g.tmax() + tau) / this->Ts);
        int nfrom = max(1, A);
        int nto = min((int)this->rmemsize,B);
        double startt = nfrom * this->Ts - tau;
        int ifrom = (int)round(( startt - this->g.tmin() ) / stepwidth);
        int istep = tabls;
        int mfrom = nfrom - 1;
        int mto = nto - 1;
        int mstep = 1;
        
        return mainFilterLoop(mfrom,mto,mstep,ifrom,istep,this->rmem,&pulsetable[0]);
        
    }

    virtual string name() const {
        return "TabulatedOerderMeyrAndMassey";
    }

protected:
    /** Table stores values of the transmit pulse g.pulse */
    std::vector<VALTYPE> pulsetable;
    /**  Width in t between elements in the table.DOUBLE IS DELIBERATE, DO NOT MODIFY. */
    double stepwidth;
    /** Multiplier for indexing the table */
    int tabls;
            
};

/**
 * Tabulated Oerder Meyer.  Store pulse in a lookup table for fast access.  Table is memory
 * aligned for extra speed.
 */
template <class Pulse>
class AlignedTabulatedOerderMeyerAndMassey : public OerderMeyerAndMassey<Pulse> {
    
public:

    AlignedTabulatedOerderMeyerAndMassey(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<VALTYPE> >& pilots,
            const Pulse& g,
            const VALTYPE T,
            const VALTYPE Ts,
            const VALTYPE taumin,
            const VALTYPE taumax,
            const unsigned int mintabs = DEFAULTTABLESIZE) : 
    OerderMeyerAndMassey<Pulse>(P,D,pilots,g,T,Ts,taumin,taumax) {
        tabls = (int) ceil(mintabs*Ts/T);
        stepwidth = Ts/tabls;
        //fill up the pulse table with the transmit pulse
        vector<VALTYPE> pulsetable;
        for (VALTYPE t = g.tmin(); t <= g.tmax(); t += stepwidth) pulsetable.push_back(g.pulse(t));
        //fill up the pulse tablstable with memory aligned transmit pulse
        tablstable.resize(pulsetable.size(), vector<VALTYPE>());
        for(int i = 0; i < pulsetable.size(); i++){
            for(int t = 0; t + i < pulsetable.size(); t += tabls) {
                tablstable[i].push_back(pulsetable[t+i]);
            }
        }
    }
    
    /** Inner product between the received signal and pulse g */
    virtual complex<VALTYPE> dotr(VALTYPE tau) const {
        //setup loop counters 
        int A = (int) ceil((this->g.tmin() + tau) / this->Ts);
        int B = (int) floor((this->g.tmax() + tau) / this->Ts);
        int nfrom = max(1, A);
        int nto = min((int)this->rmemsize,B);
        double startt = nfrom * this->Ts - tau; //DOUBLE IS DELIBERATE, DO NOT MODIFY.
        int ifrom = (int)round(( startt - this->g.tmin() ) / stepwidth);
        int mfrom = nfrom - 1;
        int mto = nto - 1;
        int mstep = 1;
        
        //return mainFilterLoop(mfrom,mto,mstep,ifrom,istep,this->rmem,&pulsetable[0]);
        return mainBankedFilterLoop(mfrom,mto,mstep,this->rmem,&tablstable[ifrom][0]);
        
    }

    virtual string name() const {
        return "AlignedTabulatedOerderMeyrAndMassey";
    }

protected:
    /** Store the pulse table in a way that encourages good memory locality */
    vector<vector<VALTYPE> > tablstable;
    /**  Width in t between elements in the table. DOUBLE IS DELIBERATE, DO NOT MODIFY.*/
    double stepwidth;
    /** Multiplier for indexing the table */
    int tabls;
};

}

#endif	/* OERDERMEYERMASSEY_H */

