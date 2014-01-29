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

#ifndef MSGTYPE
#define MSGTYPE double
//#define MSGTYPE float
#endif

typedef std::complex<MSGTYPE> complex;

class PhaseEstimator
{
 public:
  virtual ~PhaseEstimator() {};
  virtual complex complexGainEstimate() = 0;
  virtual MSGTYPE frequencyRateEstimate() = 0;
  virtual MSGTYPE frequencyOffsetEstimate() = 0;
  virtual MSGTYPE objectiveFunctionValue() = 0;
  
  //Run estimator on data y.
  virtual void estimate(const std::vector<complex>& y) = 0;
  //Allows setting data y and a new set of pilots symbols.
  virtual void estimate(const std::vector<complex>& y, const std::vector<complex>& p) = 0;

};

class CoherentMackenthun : public PhaseEstimator
{

public:
  //contruct by value
  CoherentMackenthun(
		     const std::vector<int>& D,
		     const std::vector<int>& P,
		     const std::vector<complex>& p,
		     const unsigned int M); 

  virtual void estimate(const std::vector<complex>& y);
  virtual void estimate(const std::vector<complex>& y, const std::vector<complex>& p);
  
  virtual inline complex complexGainEstimate() { return chat;}
  virtual MSGTYPE frequencyRateEstimate() { throw "No frequency rate estimate available"; }
  virtual MSGTYPE frequencyOffsetEstimate() { throw "No frequency offset estimate available"; }
  virtual MSGTYPE objectiveFunctionValue() { return Qhat; }

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
  const std::vector<complex> p; //array of pilot symbols  

  //estimator output variable
  complex chat;

  //working memory
  std::vector<IndexedReal> z;
  std::vector<std::complex<MSGTYPE> > g;

  //useful variables
  const MSGTYPE w; // = 2.0*scala.math.Pi/M
  const complex eta; // = exp(jw)
  const complex nu; // = eta - 1

 private:
  MSGTYPE Qhat; //value of objective function
  
};



class CoherentMackenthunWithDoppler : public CoherentMackenthun {

public:
  CoherentMackenthunWithDoppler(
		     const std::vector<int>& D,
		     const std::vector<int>& P,
		     const std::vector<complex>& p,
		     const unsigned int M,
		     const MSGTYPE fmin,
		     const MSGTYPE fmax,
		     const MSGTYPE T,
		     const MSGTYPE search_oversample=2.0); 
  
  //virtual void estimate(const std::vector<complex>& y);
  virtual void estimate(const std::vector<complex>& y, const std::vector<complex>& p);
  virtual inline MSGTYPE frequencyOffsetEstimate() { return fhat; }
  virtual inline MSGTYPE objectiveFunctionValue() { return Qhat; }

 protected:
  MSGTYPE       fhat; //frequency estimate
  const MSGTYPE Ts; //symbol period
  MSGTYPE fstep;
  const MSGTYPE fmin;
  const MSGTYPE fmax;
  std::vector<std::complex<MSGTYPE> > yf; //working memory

  ///Refine the Doppler estimate using Newton Raphson
  virtual void refine(const std::vector<complex>& y, const std::vector<complex>& p);
  virtual void hardDecisionsAndDerotate(const std::vector<complex>& y, const std::vector<complex>& p);
  std::vector<complex> X; //working memory for Newton Raphson

private:
  MSGTYPE Qhat; //value of objective function
  
};


class CoherentMackenthunWithDopplerAndDopplerRate : public CoherentMackenthunWithDoppler {

 public:
  CoherentMackenthunWithDopplerAndDopplerRate(
		     const std::vector<int>& D,
		     const std::vector<int>& P,
		     const std::vector<complex>& p,
		     const unsigned int M,
		     const MSGTYPE fmin,
		     const MSGTYPE fmax,
		     const MSGTYPE frmin,
		     const MSGTYPE frmax,
		     const MSGTYPE T,
		     const MSGTYPE search_oversample=2.0); 
  
  //virtual void estimate(const std::vector<complex>& y);
  virtual void estimate(const std::vector<complex>& y, const std::vector<complex>& p);
  virtual inline MSGTYPE frequencyRateEstimate() { return frhat; }
  virtual inline MSGTYPE objectiveFunctionValue() { return Qhat; }

 protected:
  MSGTYPE       frhat; //frequency rate estimate
  MSGTYPE frstep;
  const MSGTYPE frmin;
  const MSGTYPE frmax;

  ///Refine the Doppler estimate using Newton Raphson
  virtual void refine(const std::vector<complex>& y, const std::vector<complex>& p);
  virtual void hardDecisionsAndDerotate(const std::vector<complex>& y, const std::vector<complex>& p);


private:
  MSGTYPE Qhat; //value of objective function
};

#endif
