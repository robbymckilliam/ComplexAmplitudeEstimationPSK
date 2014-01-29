/* 
 * File:   TabulatedFastAlgorithm.h
 * Author: Robby McKilliam
 *
 * Created on 20 February 2013, 6:36 PM
 */

#ifndef TABULATEDFASTALGORITHM_H
#define	TABULATEDFASTALGORITHM_H

#include "TimeOffsetEstimator.h"
#include "FinitePulse.h"
#include "MainLoop.h"
#include <complex>
#include <vector>
#include <algorithm>

#ifndef VALTYPE
#define VALTYPE double
#endif
#ifndef MSGTYPE
#define MSGTYPE VALTYPE
#endif

#ifndef DEFAULTTABLESIZE
#define DEFAULTTABLESIZE 50000
#endif

namespace TimeOffset {

/**
 * Version of the polyphase TimeOffset estimator that tabulates the transmit pulse in a way to
 * allow fast access, but also fast, well pipeliped loops for filters.
 */
template <class Pulse>
class TabulatedFastAlgorithm : public FastAlgorithm<Pulse> {
    using FastAlgorithm<Pulse>::g;
    using FastAlgorithm<Pulse>::gfunc;
    using FastAlgorithm<Pulse>::a;
    using FastAlgorithm<Pulse>::b;
    using FastAlgorithm<Pulse>::n0;
    using FastAlgorithm<Pulse>::T;
    using FastAlgorithm<Pulse>::taumin;
    using FastAlgorithm<Pulse>::taumax;
    using FastAlgorithm<Pulse>::Ts;
    using FastAlgorithm<Pulse>::Delta;
    using FastAlgorithm<Pulse>::rmem;
    using FastAlgorithm<Pulse>::rmemsize;
    
public:

    TabulatedFastAlgorithm(const vector<int>& P,
		  const vector<int>& D,
		  const vector<complex<VALTYPE> >& pilots,
		  const Pulse& g,
		  const VALTYPE T,
		  const VALTYPE Ts,
		  const VALTYPE taumin,
		  const VALTYPE taumax,
		  const unsigned int c,
		  const unsigned int p,
		  const unsigned int q,
                              const unsigned int mintabs = DEFAULTTABLESIZE,
		  const VALTYPE brenttol = DEFAULTBRENTTOL) :
  FastAlgorithm<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, brenttol) {
      //int d = lcm(q,c);
      int d = c * b; //we require the number of table elements per symbol to be a multiple of d
      tabls = (mintabs+1)/d;
      int tablesizepersymbol = tabls*d;
      stepwidth = T/tablesizepersymbol;
      for(double t = g.tmin(); t <= g.tmax(); t += stepwidth) pulsetable.push_back(g.pulse(t));
    };
    
    /** Inner product between the received signal and pulse g */
    inline virtual complex<VALTYPE> dotr(VALTYPE tau) const {
        
        //step loop counters. DOUBLES ARE DELIBERATE, DO NOT MODIFY
        int A = (int) ceil((g.tmin() + tau) / Ts);
        int B = (int) floor((g.tmax() + tau) / Ts);
        int nfrom = max(1, A);
        int nto = min((int)rmemsize,B);
        double startt = nfrom * Ts - tau;
        int ifrom = (int)round(( startt - g.tmin() ) / stepwidth);
        int istep = a*tabls;
        int mfrom = nfrom - 1;
        int mto = nto - 1;
        int mstep = 1;
        
        return mainFilterLoop(mfrom, mto, mstep, ifrom, istep, rmem, &pulsetable[0]);
    }

    /** The h_{\ell,k} sequence from the paper, computed by convolution*/
    inline virtual complex<VALTYPE> h(int ell, int k) const {
         
        //step loop counters. DOUBLES ARE DELIBERATE, DO NOT MODIFY
        double A = 1 - (g.tmax() + taumin) / Delta + ((double) ell) / b;
        double B = 1 - (g.tmin() + taumin) / Delta + ((double) ell) / b;
        int Bprime = (int) ceil((k - B + ell * n0) / a);
        int Aprime = (int) floor((k - A + ell * n0) / a);
        int Bpp = a * Bprime - ell * n0;
        int App = a * Aprime - ell * n0;
        int nfrom = max((int)(ell/b+1),Bpp);
        //int nto = min((int)(a*(rmemsize+1)/b),App);
	//int nto = min((int)(a*rmemsize/b),App);
	int nto = App;
        double startt = -(k-nfrom-1)*Delta - taumin + (ell*Ts)/a; 
        int ifrom = (int)round(( startt - g.tmin() ) / stepwidth);
        int istep = a*b*tabls;
        int mfrom = (b*nfrom+ell)/a-1;
        //int mto = (b*nto+ell)/a-1;
        int mto = min((int)rmemsize,(int)((b*nto+ell)/a))-1;
	int mstep = b;
        
        return mainFilterLoop(mfrom, mto, mstep, ifrom, istep, rmem, &pulsetable[0]);
    }
    
protected:
    /** Table stores values of the transmit pulse g.pulse */
    vector<VALTYPE> pulsetable;
    /** 
     * Width in t between elements in the table.
     * THIS IS DELIBERATELY A DOUBLE FOR PRECISION. DO NOT CHANGE THIS!!!
     */
    double stepwidth;
    /** Multiplier for indexing the table */
    int tabls;
    
};

template <class Pulse>
class TabulatedDirect : public TabulatedFastAlgorithm<Pulse> {
    
public:

    TabulatedDirect(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<VALTYPE> >& pilots,
            const Pulse& g,
            const VALTYPE T,
            const VALTYPE Ts,
            const VALTYPE taumin,
            const VALTYPE taumax,
            const unsigned int c,
            const unsigned int p,
            const unsigned int q,
            const unsigned int mintabs = DEFAULTTABLESIZE,
            const VALTYPE brenttol = DEFAULTBRENTTOL) :
  TabulatedFastAlgorithm<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, mintabs, brenttol) {
    }

    /** Return the values of Z computed on the grid taumin to taumax by Ts/c.  Direct computation. */
    virtual void fillZgrid() {
        for (unsigned int k = 1; k <= this->K; k++) {
            VALTYPE sum = 0.0;
            for(unsigned int i = 0; i < this->D.size(); i++) sum += std::abs(this->bfunc(k + this->c * (this->D[i] + 1)));
            //for (int i : this->D) sum += std::abs(this->bfunc(k + this->c * (i + 1)));
            this->Zgrid[k - 1] = sum;
        }
    }

    /* Return the values of Y computed on the grid taumin to taumax by T/c. Direct convolution */
    virtual void fillYgrid() {
        for (unsigned int k = 1; k <= this->K; k++) {
            complex<VALTYPE> sum = 0.0;
            for (unsigned int i = 0; i < this->P.size(); i++)
                sum += this->bfunc(k + this->c * (this->P[i] + 1)) * conj(this->pilots[i]);
            this->Ygridmag[k - 1] = std::abs(sum);
        }
    }

    virtual string name() const {
        return "TabulatedDirect";
    }
};


/**
 * Recursive estimator of time offset.  The data symbol indices must now be 
 * contiguous.  Exception is thrown if they are not.
 */
template <class Pulse>
class TabulatedRecursive : public TabulatedDirect<Pulse> {
public:
    
    TabulatedRecursive(const vector<int>& P,
	      const vector<int>& D,
	      const vector<complex<VALTYPE> >& pilots,
	      const Pulse& g,
	      const VALTYPE T,
	      const VALTYPE Ts,
	      const VALTYPE taumin,
	      const VALTYPE taumax,
	      const unsigned int c,
	      const unsigned int p,
	      const unsigned int q,
                    const unsigned int mintabs = DEFAULTTABLESIZE,
	      const VALTYPE brenttol = DEFAULTBRENTTOL) :
  TabulatedDirect<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q,mintabs,brenttol) {
        //check D is contiguous and sorted
        for (unsigned int i = 0; i < D.size() - 1; i++) {
            if (D[i] != (D[i + 1] - 1)) {
                std::ostringstream str;
                str << "Exception constructing Recursive. The data symbols must be contiguous." << endl;
                throw str.str();
            }
        }
    }

    /* 
     * Return the values of Z computed on the grid taumin to taumax by T/c. 
     * Uses recursive algorithm that only applies when D is contiguous.
     */
    virtual void fillZgrid() {
        for (unsigned int k = 1; k <= this->c; k++) {
            VALTYPE sum = 0.0;
            for (unsigned int i = 0; i < this->D.size(); i++) sum += abs(this->bfunc(k + this->c * (this->D[i] + 1)));
            this->Zgrid[k - 1] = sum;
            for (unsigned int m = 0; k + (m + 1) * this->c <= this->K; m++)
                this->Zgrid[k - 1 + (m + 1) * this->c] = this->Zgrid[k - 1 + m * this->c] - abs(this->bfunc(k + m * this->c + (this->Dmin + 1) * this->c)) + abs(this->bfunc(k + m * this->c + (this->Dmax + 1) * this->c + this->c));
        }
    }

    virtual string name() const {
        return "TabulatedRecursive";
    }
};

}

#endif	/* TABULATEDFASTALGORITHM_H */

