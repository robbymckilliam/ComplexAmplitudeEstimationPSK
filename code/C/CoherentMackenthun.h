//
//  CoherentMackenthun.h
//  CoherentMackenthun
//
//  Created by Robby McKilliam on 4/06/12.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#ifndef _CoherentMackenthun_h_included_
#define _CoherentMackenthun_h_included_

#include "IndexedReal.h"
#include <complex>
#include <vector>

typedef std::complex<double> complexd;

class PhaseEstimator
{
 public:
  virtual ~PhaseEstimator() {};
  virtual complexd complexGainEstimate() = 0;
  virtual double objectiveFunctionValue() = 0;
  virtual double noiseVarianceEstimate() = 0;
  
  //Run estimator on data y.
  virtual void estimate(const std::vector<complexd>& y) = 0;
  //Allows setting data y and a new set of pilots symbols.
  virtual void estimate(const std::vector<complexd>& y, const std::vector<complexd>& p) = 0;

};

class CoherentMackenthun : public PhaseEstimator
{

public:
  //construct by value
  CoherentMackenthun(
		     const std::vector<int>& D,
		     const std::vector<int>& P,
		     const std::vector<complexd>& p,
		     const unsigned int M); 

  virtual void estimate(const std::vector<complexd>& y);
  virtual void estimate(const std::vector<complexd>& y, const std::vector<complexd>& p);
  
  virtual inline complexd complexGainEstimate() { return chat;}
  virtual double noiseVarianceEstimate() { return noisevarhat; }
  virtual double objectiveFunctionValue() { return Qhat; }

  const unsigned int L; //total number of pilots
  const unsigned int M; //M for M-PSK
  const unsigned int absD; //number of data symbols
  const unsigned int absP; //number of data symbols
    
protected:
  
  inline int sigma(int k) { return z[k % absD].i; }
  
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

class PerfectChannel : public PhaseEstimator {
    
public:
    const complexd hatc;
    const double variance;
    
    PerfectChannel(const complexd c, const double v) : hatc(c), variance(v) {}
    
    virtual void estimate(const std::vector<complexd>& y) {}
    virtual void estimate(const std::vector<complexd>& y, const std::vector<complexd>& p) {}
    virtual inline complexd complexGainEstimate() { return hatc;}
    virtual double noiseVarianceEstimate() { return variance; }
    virtual double objectiveFunctionValue() { return 0.0; }
    
};

#endif
