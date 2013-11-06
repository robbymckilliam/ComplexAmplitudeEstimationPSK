/* 
 * File:   FinitePulse.h
 * Author: Robby McKilliam
 *
 * Created on 23 January 2013, 1:10 PM
 */

#ifndef FINITEPULSE_H
#define	FINITEPULSE_H

#include <vector>
#include <iostream>
#include "Util.h"

#ifndef VALTYPE
#define VALTYPE double
#endif
#ifndef MSGTYPE
#define MSGTYPE VALTYPE
#endif

using namespace std;

namespace TimeOffset {

class FinitePulse {
public:
    virtual ~FinitePulse() {
    };
    /** pulse is zero for t less than tmin */
    virtual VALTYPE tmin() const = 0;
    /** pulse is zero for t more than tmax */
    virtual VALTYPE tmax() const = 0;
    /** pulse period */
    virtual VALTYPE T() const = 0;
    /** The transmit pulse */
    virtual VALTYPE pulse(VALTYPE t) const = 0;
};

/** Truncated sinc pulse goes out to the numzeros zero in both positive and negative directions */
class TruncatedSincPulse : public FinitePulse {
public:
    
    TruncatedSincPulse(VALTYPE T, unsigned int numzeros) : tmin_(-T*numzeros), tmax_(T*numzeros), T_(T), sqrtT(sqrt(T)) {};
 TruncatedSincPulse(VALTYPE T, VALTYPE duration, bool dodgyflag) : tmin_(-T*(int)ceil(duration/T/2)), tmax_(T*(int)ceil(duration/T/2)), T_(T), sqrtT(sqrt(T)) {};
    //TruncatedSincPulse() = delete; //no default constructor //c++11 only

    /** The truncated sinc pulse */
    inline VALTYPE pulse(VALTYPE t) const {
        if (t > tmin_ && t < tmax_) return sinc(t / T_) / sqrtT;
        else return 0.0;
    }

    inline VALTYPE tmin() const {
        return tmin_;
    }

    inline VALTYPE tmax() const {
        return tmax_;
    }

    inline VALTYPE T() const {
        return T_;
    }

protected:
    const VALTYPE tmin_;
    const VALTYPE tmax_;
    const VALTYPE T_;
    const VALTYPE sqrtT;
private:

    TruncatedSincPulse() : tmin_(0), tmax_(0), T_(0), sqrtT(0) {
    }; //use delete instead for c++11
};

/** 
 * Takes a pulse of finite duration and normalises it to have energy 1.
 * Uses numerical integration to compute the normalising constant
 */
template <class P> //P should be a FinitePulse
class NormalisedFinitePulse : public FinitePulse {
public:

    NormalisedFinitePulse(const P& p_) : p(p_), tmin_(p_.tmin()), tmax_(p_.tmax()), T_(p_.T()) {
        VALTYPE energy = trapezoidal(normpulse, tmin_, tmax_, 100000); //last number is integration steps
        normalisingconstant = sqrt(energy);
    };
    //~NormalisedFinitePulse() { delete p; }
    //NormalisedFinitePulse( const NormalisedFinitePulse& other ) = delete; //no copy
    //NormalisedFinitePulse& operator =(const NormalisedFinitePulse&) = delete; //no assignment
    //NormalisedFinitePulse() = delete; //no default constructor

    VALTYPE normpulse(VALTYPE t) {
        VALTYPE d = p.pulse(t);
        return d*d;
    }

    inline VALTYPE tmin() const {
        return tmin_;
    }

    inline VALTYPE tmax() const {
        return tmax_;
    }

    inline VALTYPE T() const {
        return T_;
    }

    /** The normalised pulse */
    inline VALTYPE pulse(VALTYPE t) const {
        return p.pulse(t) / normalisingconstant;
    }

protected:
    const P p;
    const VALTYPE tmin_;
    const VALTYPE tmax_;
    const VALTYPE T_;
    VALTYPE normalisingconstant;
    
private:
    NormalisedFinitePulse() : p(0), tmin_(0), tmax_(0), T_(0) {};//use delete in c++11

};

class TruncatedRootRaisedCosine : public FinitePulse {
public:

    TruncatedRootRaisedCosine(VALTYPE T, VALTYPE beta, VALTYPE duration, bool dodgyflag) : T_(T), beta(beta)  {
        const VALTYPE stepsize = 0.01; //step taken whilst looking for zeros (this will work only if zeros are atleast stepsize apart)
        VALTYPE c = 0.0;
        int csign = 1;
        while(c < duration/T/2) {
            while (signum(rootraisedcosine(c)) == csign) c += stepsize;
            csign = -csign;
        }
        RealRealMemberFunctionCaller<TruncatedRootRaisedCosine> f(this, &TruncatedRootRaisedCosine::rootraisedcosine);
        VALTYPE nt = Bisection<RealRealMemberFunctionCaller<TruncatedRootRaisedCosine> >(f, c - stepsize, c, 1e-6, 100).zero();
        tmin_ = -nt*T;
        tmax_ = nt*T;
        VALTYPE dur = tmax_ - tmin_;
        //std::cout << dur << " but request duration is " << duration << std::endl;
        if( dur > 2*duration || dur < duration/2 ) { //check that the duration obtained is reasonable
            std::ostringstream str;
            str << "Something when wrong with the bisection method, duration is ";
            str << dur << " but request duration is " << duration;
            throw str.str();
        }
    }
    TruncatedRootRaisedCosine(VALTYPE T, VALTYPE beta, int numzeros) : T_(T), beta(beta) {
        const VALTYPE stepsize = 0.01; //step taken whilst looking for zeros (this will work only if zeros are atleast stepsize apart)
        VALTYPE c = 0.0;
        int csign = 1;
        for (int i = 1; i <= numzeros; i++) {
            while (signum(rootraisedcosine(c)) == csign) c += stepsize;
            csign = -csign;
        }
        RealRealMemberFunctionCaller<TruncatedRootRaisedCosine> f(this, &TruncatedRootRaisedCosine::rootraisedcosine);
        VALTYPE nt = Bisection<RealRealMemberFunctionCaller<TruncatedRootRaisedCosine> >(f, c - stepsize, c, 1e-6, 100).zero();
        tmin_ = -nt*T;
        tmax_ = nt*T;
    }
    //TruncatedRootRaisedCosine() = delete; //no default constructor //c++11 only

    inline VALTYPE pulse(VALTYPE t) const {
        if (t > tmin_ && t < tmax_) return rootraisedcosine(t / T_);
        else return 0.0;
    }

    inline VALTYPE tmin() const {
        return tmin_;
    }

    inline VALTYPE tmax() const {
        return tmax_;
    }

    inline VALTYPE T() const {
        return T_;
    }

    /** A root raised cosine with period 1 and rolloff beta */
    VALTYPE rootraisedcosine(VALTYPE t) const {
        VALTYPE abst = fabs(t); //pulse is symmetric about zero, so just use magnitude of t
        if (abst < 5e-3) { //second order expansion if t is near zero
            VALTYPE term0 = 1 + beta * (4 / pi - 1);
            VALTYPE term2 = (cub((beta - 1) * pi) + 96 * beta * beta * (4 * beta + pi - beta * pi) - 12 * beta * sqr(pi + beta * pi)) / 6 / pi;
            return term0 + term2 * abst*abst;
        }
        if (fabs(abst - 1.0 / 4 / beta) < 5e-4) { //first order expansion if t is near 1/(4beta)
            VALTYPE a = (1 + beta) * pi / 4 / beta;
            VALTYPE term0 = beta * (sin(a) - 2 * cos(a) / pi);
            VALTYPE term1 = beta * ((12 * beta + pi * pi) * cos(a) + 2 * (1 - 2 * beta) * pi * sin(a)) / pi;
            return term0 + term1 * (abst - 1.0 / 4 / beta);
        } else { //otherwise use direct formula
            VALTYPE a = pi*t;
            VALTYPE b = beta*t;
            return (sin(a * (1 - beta)) + 4 * b * cos(a * (1 + beta))) / (a * (1 - 16 * b * b));
        }
    }

protected:
    VALTYPE tmin_;
    VALTYPE tmax_;
    const VALTYPE T_;
    const VALTYPE beta;
private:

    TruncatedRootRaisedCosine() : tmin_(0), tmax_(0), T_(0), beta(0) {
    }; //use delete in c++11
};

/** Store the values of the pulse in a lookup table for speed */
template <class P> //P should be a FinitePulse
class FinitePulseWithLookupTable : public FinitePulse {
public:

    FinitePulseWithLookupTable(const P& p_, const unsigned int tablesize_) :
    p(p_), T_(p_.T()), tmin_(p_.tmin()), tmax_(p_.tmax()), tablesize(tablesize_), stepwidth((tmax_ - tmin_) / (tablesize_ - 1)) {
        for (int i = 0; i < tablesize; i++) table.push_back(p.pulse(tmin_ + i * stepwidth));
    }

    virtual VALTYPE tmin() const {
        return tmin_;
    }

    virtual VALTYPE tmax() const {
        return tmax_;
    }

    virtual VALTYPE T() const {
        return T_;
    }

    virtual VALTYPE pulse(VALTYPE t) const {
        if (t < tmax_ && t > tmin_) {
            int ai = (int) floor((t - tmin_) / stepwidth);
            int bi = (int) ceil((t - tmin_) / stepwidth);
            VALTYPE x = t - tmin_ - stepwidth*ai;
            VALTYPE m = (table[bi] - table[ai]) / stepwidth;
            return m * x + table[ai];
        } else return 0.0;
    }

protected:
    const P p;
    const unsigned int tablesize;
    vector<VALTYPE> table;
    const VALTYPE tmin_;
    const VALTYPE tmax_;
    const VALTYPE T_;
    const VALTYPE stepwidth;

};

}

#endif	/* FINITEPULSE_H */

