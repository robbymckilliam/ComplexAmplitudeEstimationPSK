/* 
 * File:   FastAlgorithm.h
 * Author: Robby McKilliam
 *
 * Created on 4 February 2013, 8:13 PM
 */

#ifndef FASTALGORITHM_H
#define	FASTALGORITHM_H

#include "TimeOffsetEstimator.h"
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

namespace TimeOffset {


/**
 * Override Naive with a faster way to compute the b vectors.  This does not store the bk sequence.
 */
template <class Pulse> //P should be a FinitePulse
class FastAlgorithmNoStoreb : public Naive<Pulse> {
    
public:

    FastAlgorithmNoStoreb(const vector<int>& P,
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
			  const VALTYPE brenttol = DEFAULTBRENTTOL) :
  Naive<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, brenttol),
    Dmin(*std::min_element(D.begin(), D.end())),
    Dmax(*std::max_element(D.begin(), D.end())),
    Pmin(*std::min_element(P.begin(), P.end())),
    Pmax(*std::max_element(P.begin(), P.end())),
    p(p),
    q(q),
    K(floor((taumax - taumin) / this->Delta)),
    Zgrid(K, 0.0),
    Ygridmag(K, 0.0),
    SSgrid(K, 0.0) {
        //setup a, b, n0 and m0 for the polyphase filter
        int d = gcd(q, c * p);
        a = c * p / d;
        b = q / d;
        extended_gcd(b, a, d, n0, m0);
    }

    virtual VALTYPE coarseMaximiseSS() {
        fillSSgrid(); //compute objective function on the grid
        std::vector<VALTYPE>::iterator biggest = std::max_element(SSgrid.begin(), SSgrid.end());
        int khat = std::distance(SSgrid.begin(), biggest); //get position of max element
        //std::cout << "dir or rec " <<  khat << " and " << SSgrid[khat] << std::endl;
        return this->taumin + khat*this->Delta;
    }

    /** fills the vector SSgrid with values of the objective function SS */
    void fillSSgrid() {
        fillZgrid();
        fillYgrid();
        for (unsigned int i = 0; i < K; i++) SSgrid[i] = Zgrid[i] + Ygridmag[i];
    }

    /** zero filled received signal z_n in the paper */
    inline complex<VALTYPE> zn(int n) const {
        if (n % a == 0) return this->r(n / a);
        else return complex<VALTYPE>(0, 0);
    }

    /** banked received signal z_{ell,n} from the paper*/
    inline complex<VALTYPE> z(int ell, int n) const {
        return zn(b * n + ell);
    }

    /** The sequence g_{ell,n} from the paper */
    virtual VALTYPE gfunc(int ell, int n) const {
        //return conj(this->g.pulse(-(n - 1) * this->Delta - this->taumin + (ell * this->Ts) / a));
        return this->g.pulse(-(n - 1) * this->Delta - this->taumin + (ell * this->Ts) / a); //pulse assumed real conjugate not necessary
    }

    /** The sequence bk from the paper */
    virtual complex<VALTYPE> bfunc(int k) const {
        complex<VALTYPE> sum(0, 0);
        for (unsigned int ell = 0; ell < b; ell++) sum += h(ell, k);
        //cout << b << ", " << sum << endl;
        return sum;
    }

    /** The h_{\ell,k} sequence from the paper, computed by convolution*/
    virtual complex<VALTYPE> h(int ell, int k) const {
        VALTYPE A = 1 - (this->g.tmax() + this->taumin) / this->Delta + ((VALTYPE) ell) / b;
        VALTYPE B = 1 - (this->g.tmin() + this->taumin) / this->Delta + ((VALTYPE) ell) / b;
        int Bprime = (int) ceil((k - B + ell * n0) / a);
        int Aprime = (int) floor((k - A + ell * n0) / a);
        complex<VALTYPE> sum(0, 0);
        int Bpp = a * Bprime - ell * n0;
        int App = a * Aprime - ell * n0;
        for (int n = Bpp; n <= App; n += a) sum += z(ell, n) * gfunc(ell, k - n);
        return sum;
    }

    /** fill vector Zgrid with discretized Z function */
    virtual void fillZgrid() = 0;

    /** fill vector Ygridmag with discretized magnitude of the Y function */
    virtual void fillYgrid() = 0;

    //output for testing

    vector<VALTYPE> getZgrid() {
        return Zgrid;
    }

    vector<VALTYPE> getYgrid() {
        return Ygridmag;
    }

    vector<VALTYPE> getSSgrid() {
        return SSgrid;
    }

protected:
    /** Minimum data symbol index */
    const int Dmin;
    /** Maximum data symbol index */
    const int Dmax;
    /** Minimum pilot symbol index */
    const int Pmin;
    /** Maximum data symbol index */
    const int Pmax;
    const unsigned int p;
    const unsigned int q;
    unsigned int a;
    unsigned int b;
    int n0;
    int m0;
    /** Size of the search grid */
    const unsigned int K;
    vector<VALTYPE> Zgrid;
    vector<VALTYPE> Ygridmag;
    vector<VALTYPE> SSgrid;

};

/**
 * Fast polyphase computer for bk sequence.  Stores the bk sequence for fast access.
 */
template <class Pulse>
class FastAlgorithm : public FastAlgorithmNoStoreb<Pulse> {
public:

    FastAlgorithm(const vector<int>& P,
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
		  const VALTYPE brenttol = DEFAULTBRENTTOL) :
  FastAlgorithmNoStoreb<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, brenttol),
    bstore(this->K + c*(max(this->Dmax, this->Pmax) - min(this->Dmin, this->Pmin)), complex<VALTYPE>(0, 0)),
    boffset(1 + c*(min(this->Dmin, this->Pmin) + 1)) {
    };

    virtual VALTYPE estimate(const vector<complex<VALTYPE> >& r) {
        this->setupr(r);
        for (unsigned int k = 0; k < bstore.size(); k++) bstore[k] = FastAlgorithmNoStoreb<Pulse>::bfunc(k + boffset);
        VALTYPE tautilde = this->coarseMaximiseSS();
        VALTYPE tauhat = this->refineCoarseEstimate(tautilde);
        return tauhat;
    }
    
    //b takes from array
    inline virtual complex<VALTYPE> bfunc(int k) const {
        return bstore[k - boffset];
    }

protected:
    vector<complex<VALTYPE> > bstore;
    const int boffset;

};

template <class Pulse>
class Direct : public FastAlgorithm<Pulse> {
public:

    Direct(const vector<int>& P,
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
            const VALTYPE brenttol = DEFAULTBRENTTOL) :
  FastAlgorithm<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, brenttol) {
    }

    /** Return the values of Z computed on the grid taumin to taumax by Ts/c.  Direct computation. */
    virtual void fillZgrid() {
        for (unsigned int k = 1; k <= this->K; k++) {
            VALTYPE sum = 0.0;
            for(int i = 0; i < this->D.size(); i++) sum += std::abs(this->bfunc(k + this->c * (this->D[i] + 1)));
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
        return "Direct";
    }
};


/**
 * Recursive estimator of time offset.  The data symbol incdices must now be 
 * contiguous.  Exception is thrown if they are not.
 */
template <class Pulse>
class Recursive : public Direct<Pulse> {
public:
    
    Recursive(const vector<int>& P,
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
	      const VALTYPE brenttol = DEFAULTBRENTTOL) :
  Direct<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, brenttol) {
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
            for (int i = 0; i < this->D.size(); i++) sum += abs(this->bfunc(k + this->c * (this->D[i] + 1)));
            this->Zgrid[k - 1] = sum;
            for (int m = 0; k + (m + 1) * this->c <= this->K; m++)
                this->Zgrid[k - 1 + (m + 1) * this->c] = this->Zgrid[k - 1 + m * this->c] - abs(this->bfunc(k + m * this->c + (this->Dmin + 1) * this->c)) + abs(this->bfunc(k + m * this->c + (this->Dmax + 1) * this->c + this->c));
        }
    }

    virtual string name() const {
        return "Recursive";
    }
};

}

#endif	/* FASTALGORITHM_H */

