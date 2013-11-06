/* 
 * File:   TimeOffsetEstimator.h
 * Author: Robby McKilliam
 *
 * Created on 23 January 2013, 10:58 AM
 */

#ifndef TIMEOFFSETESTIMATOR_H
#define	TIMEOFFSETESTIMATOR_H

#include "FinitePulse.h"
#include <complex>
#include <vector>
#include <algorithm>

#ifndef VALTYPE
#define VALTYPE double
#endif
#ifndef MSGTYPE
#define MSGTYPE VALTYPE
#endif

#define DEFAULTBRENTTOL 1e-6

using namespace std;

namespace TimeOffset {

/** Interface for time offset estimators  */
class TimeOffsetEstimator {
public:

    virtual ~TimeOffsetEstimator() {
    };
    /** Run the estimator. Return time offset estimate */
    virtual VALTYPE estimate(const vector<complex<VALTYPE> >& r) = 0;
    /** Return a string with the name of this estimator */
    virtual string name() const = 0;
};

/** Naive implementation of the time offset estimator that computes the objective function directly */
template <class Pulse> //P should be a FinitePulse
class Naive : public TimeOffsetEstimator {
public:

    Naive(const vector<int>& P_,
	  const vector<int>& D_,
	  const vector<complex<VALTYPE> >& pilots_,
	  const Pulse& g_,
	  const VALTYPE T_,
	  const VALTYPE Ts_,
	  const VALTYPE taumin_,
	  const VALTYPE taumax_,
	  const unsigned int c_,
	  const VALTYPE brenttol_ = DEFAULTBRENTTOL) :
  L(D_.size() + P_.size()),
    Delta(T_ / c_),
    P(P_),
    D(D_),
    pilots(pilots_),
    g(g_),
    T(T_),
    Ts(Ts_),
    taumin(taumin_),
    taumax(taumax_),
    c(c_),
    brenttol(brenttol_*T) { 
        if(taumax <= taumin) throw "Maximum time offset taumax must be larger than minimum time offset taumin";
  };

    virtual VALTYPE estimate(const vector<complex<VALTYPE> >& r) {
        setupr(r);
        VALTYPE tautilde = coarseMaximiseSS();
        VALTYPE tauhat = refineCoarseEstimate(tautilde);
        return tauhat;
    }
    
    virtual void setupr(const vector<complex<VALTYPE> >& r) {
        rmem = r.data();
        rmemsize = r.size();
    }

    /** Obtain a coarse estimate of the time offset */
    virtual VALTYPE coarseMaximiseSS() {
        VALTYPE maxSS = -1.0; //SS will always be positive, so this is fine.
        VALTYPE taubest = taumin;
        for (VALTYPE tau = taumin; tau <= taumax; tau += Delta) {
            VALTYPE thisSS = SS(tau);
            if (maxSS < thisSS) {
                maxSS = thisSS;
                taubest = tau;
            }
        }
        //std::cout << SS(taubest) << std::endl;
        return taubest;
    }

    /** Refines the coarse estimate. */
    VALTYPE refineCoarseEstimate(VALTYPE tautilde) const {
        VALTYPE a = tautilde - Delta;
        VALTYPE c = tautilde + Delta;
        //std::cout << a - c << std::endl;
        //auto f = [this] (VALTYPE tau){return -this->SS(tau);}; //function to minimise //c++11 only
        RealRealMemberFunctionCaller<Naive> f(this,&Naive::negativeSS);
        Brent<RealRealMemberFunctionCaller<Naive> > opt(f, a, tautilde, c, brenttol);
        return opt.xmin();
    }
    VALTYPE negativeSS(VALTYPE tau) const { return -SS(tau); }

    /** Inner product between the received signal and pulse g */
    virtual complex<VALTYPE> dotr(VALTYPE tau) const {
        int A = (int) ceil((g.tmin() + tau) / Ts);
        int B = (int) floor((g.tmax() + tau) / Ts);
        complex<VALTYPE> sum(0, 0);
        //for (int n = A; n <= B; n++) sum += r(n) * conj(g.pulse(n * Ts - tau));
        for (int n = A; n <= B; n++) sum += r(n) * g.pulse(n * Ts - tau); //pulse assumed real conjugate not necessary
        return sum;
    }

    /** The Y function computing a correlation with the data and pilots */
    complex<VALTYPE> Y(VALTYPE tau) const {
        complex<VALTYPE> sum(0, 0);
        for (unsigned int i = 0; i < P.size(); i++) sum += dotr((P[i] + 1) * T + tau) * conj(pilots[i]);
        return sum;
    }

    /** The amplitude accumulating Z function */
    VALTYPE Z(VALTYPE tau) const {
        VALTYPE sum = 0.0;
        for (unsigned int i = 0; i < D.size(); i++) sum += abs(dotr((D[i] + 1) * T + tau));
        //for (int i : D) sum += abs(dotr((i + 1) * T + tau)); //c++11 only
        return sum;
    }

    /** The objective function */
    inline VALTYPE SS(VALTYPE tau) const {
        return Z(tau) + abs(Y(tau));
    }

    virtual string name() const {
        return "Naive";
    }

protected:

    /** Pointer to current received signal */
    const complex<VALTYPE>* rmem;
    /** size of current recieved signal */
    unsigned int rmemsize;
    /** Total number of transmitted symbols */
    const unsigned int L;
    /** Grid search width */
    const VALTYPE Delta;

    /** The received sequence, zeros returned for unknown values */
    inline complex<VALTYPE> r(int n) const {
      if (n > 0 && n <= (int)rmemsize) return rmem[n - 1];
        else return complex<VALTYPE>(0, 0);
    }

    const vector<int> P;
    const vector<int> D;
    const vector<complex<VALTYPE> > pilots;
    const Pulse g;
    const VALTYPE T;
    const VALTYPE Ts;
    const VALTYPE taumin;
    const VALTYPE taumax;
    const unsigned int c;
    const VALTYPE brenttol;
};

}

#endif	/* TIMEOFFSETESTIMATOR_H */

