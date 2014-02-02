//
//  CoherentMackenthun.h
//  CoherentMackenthun
//
//  Created by Robby McKilliam on 4/06/12.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#ifndef _CoherentMackenthun_h_included_
#define _CoherentMackenthun_h_included_

#include "Util.h"
#include <complex>
#include <vector>
#include <algorithm>

/** Base class for phase estimators */
class PhaseEstimator {
public:

    virtual ~PhaseEstimator() {
    };
    virtual complexd complexGainEstimate() const = 0;
    virtual double objectiveFunctionValue() const = 0;
    ///return an estimate of the variance of the complex noise
    virtual double noiseVarianceEstimate() const = 0;

    //Run estimator on data y.
    virtual void estimate(const std::vector<complexd>& y) = 0;
    //Allows setting data y and a new set of pilots symbols.
    virtual void estimate(const std::vector<complexd>& y, const std::vector<complexd>& p) = 0;

};

/** Class for sorting and storing permuations */
class IndexedReal
{

public:
  IndexedReal( double value, int index ) : v(value), i(index) {}

  bool operator<(const IndexedReal& other) const {
  return v < other.v;
  }
  
  double v;
  int i;

};

/** Semi-blind least squares estimator.  Requires O(LlogL) operations. */
class CoherentMackenthun : public PhaseEstimator {
public:
    //construct by value

    CoherentMackenthun(
            const std::vector<int>& Din,
            const std::vector<int>& Pin,
            const std::vector<complexd>& pin,
            const unsigned int Min) :
    L(Din.size() + Pin.size()),
    M(Min),
    absD(Din.size()),
    absP(Pin.size()),
    D(Din),
    P(Pin),
    p(pin),
    w(2 * pi / M),
    eta(std::polar<double>(1.0, w)),
    nu(std::polar<double>(1.0, w) - complexd(1.0, 0.0)) {
        if (pin.size() != P.size()) throw std::string("The number of pilot symbols does not match the number of indices in Pin");
        if (M < 1) throw std::string("M must be a positive integer");
        //allocate memory
        z.resize(D.size(), IndexedReal(0, 0));
        for (unsigned int i = 0; i < z.size(); i++) z[i] = IndexedReal(0, i); //initialise z with the identity permutation 
        g.resize(D.size(), complexd(0, 0));
        PUD.resize(L, 0); //union P and D
        for (unsigned int i = 0; i < P.size(); i++) PUD[i] = P[i];
        for (unsigned int i = 0; i < D.size(); i++) PUD[i + P.size()] = D[i];
        max_index = *std::max_element(PUD.begin(), PUD.end());
        min_index = *std::min_element(PUD.begin(), PUD.end());
    }

    virtual void estimate(const std::vector<complexd>& y) {
        estimate(y, p);
    }

    virtual void estimate(const std::vector<complexd>& y, const std::vector<complexd>& p) {
        if (y.size() <= max_index) {
            std::ostringstream str;
            str << "Exception in function CoherentMackenthun::estimate. The data length y is ";
            str << y.size() << " but the maximum symbol index is " << max_index;
            throw str.str();
        }
        if (p.size() != P.size()) {
            std::ostringstream str;
            str << "Exception in function CoherentMackenthun::estimate. The number of pilot symbols given is ";
            str << p.size() << " but the number of indices is " << P.size();
            throw str.str();
        }

        //the value A, norm of received signal
        double A = 0.0;
        for (int i = 0; i < L; i++) A += std::norm(y[PUD[i]]);

        //setup sequences and sort
        complexd Y(0.0, 0.0);
        for (unsigned int i = 0; i < P.size(); i++) Y += y[P[i]] * conj(p[i]);
        for (unsigned int i = 0; i < D.size(); i++) {
            double phi = std::arg(y[D[i]]);
            double u = w * round(phi / w);
            g[i] = y[D[i]] * std::polar<double>(1, -u);
            z[i] = IndexedReal(phi - u, i);
            Y += g[i];
        }
        std::sort(z.begin(), z.end());

        chat = Y / ((double) L);
        Qhat = std::norm(Y) / L;
        noisevarhat = (A - Qhat) / L;
        for (unsigned int k = 0; k < M * D.size(); k++) {
            Y += nu * g[sigma(k)];
            g[sigma(k)] *= eta;
            double Q = std::norm(Y) / L;
            if (Q > Qhat) {
                chat = Y / ((double) L);
                Qhat = Q;
                noisevarhat = (A - Qhat) / L;
            }
        }
    }

    virtual inline complexd complexGainEstimate() const {
        return chat;
    }

    virtual double noiseVarianceEstimate() const {
        return noisevarhat;
    }

    virtual double objectiveFunctionValue() const {
        return Qhat;
    }

    const unsigned int L; //total number of pilots
    const unsigned int M; //M for M-PSK
    const unsigned int absD; //number of data symbols
    const unsigned int absP; //number of data symbols

protected:

    inline int sigma(int k) {
        return z[k % absD].i;
    }

    //pointers
    const std::vector<int> D; //data positions
    const std::vector<int> P; //pilot positions
    std::vector<int> PUD; //all symbols positions (union of P and D)
    unsigned int max_index; //large symbols index, equivalent to max(PUD).
    unsigned int min_index; //large symbols index, equivalent to max(PUD).
    const std::vector<complexd> p; //array of pilot symbols  

    //estimator output variable
    complexd chat;

    //estimate of noise variance
    double noisevarhat;

    //working memory
    std::vector<IndexedReal> z;
    std::vector<std::complex<double> > g;

    //useful variables
    const double w; // = 2.0*scala.math.Pi/M
    const complexd eta; // = exp(jw)
    const complexd nu; // = eta - 1

private:
    double Qhat; //value of objective function

};

#endif
