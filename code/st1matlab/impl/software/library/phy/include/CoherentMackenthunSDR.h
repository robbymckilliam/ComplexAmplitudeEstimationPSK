//
//  CoherentMackenthun.h
//  CoherentMackenthun
//
//  Created by Robby McKilliam on 4/06/12.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#ifndef _CoherentMackenthunSDR_h_included_
#define _CoherentMackenthunSDR_h_included_

#include "IndexedRealSDR.h"
#include <complex>
#include <vector>

#ifndef CMACK_MSGTYPE
#define CMACK_MSGTYPE double
//#define CMACK_MSGTYPE float
#endif

typedef std::complex<CMACK_MSGTYPE> mycomplex;

class PhaseEstimator
{
 public:
  virtual ~PhaseEstimator() {}
  virtual mycomplex complexGainEstimate() = 0;
  virtual CMACK_MSGTYPE frequencyRateEstimate() = 0;
  virtual CMACK_MSGTYPE frequencyOffsetEstimate() = 0;
  virtual CMACK_MSGTYPE objectiveFunctionValue() = 0;
  
  //Run estimator on data y.  Will delete y when finished
  virtual void estimate(const std::vector<mycomplex>* y) = 0;
  //Run estimator on data y.
  virtual void estimate(const std::vector<mycomplex>& y) = 0;
  
  //Allows setting data y and a new set of pilots symbols.  Will delete y and p when finished
  virtual void estimate(const std::vector<mycomplex>* y, const std::vector<mycomplex>* p) = 0;
  //Allows setting data y and a new set of pilots symbols.
  virtual void estimate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p) = 0;

};

class CoherentMackenthun : public PhaseEstimator
{

public:
  //contruct by value
  CoherentMackenthun(
		     const std::vector<int>* D,
		     const std::vector<int>* P,
		     const std::vector<mycomplex>* p,
		     const int M); 
  
  ~CoherentMackenthun(void);

  virtual void estimate(const std::vector<mycomplex>* y); //will delete y when finished
  virtual void estimate(const std::vector<mycomplex>& y);
  virtual void estimate(const std::vector<mycomplex>* y, const std::vector<mycomplex>* p); //will delete y and p when finished
  virtual void estimate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p);
  
  virtual inline mycomplex complexGainEstimate() { return chat;}
  virtual CMACK_MSGTYPE frequencyRateEstimate() { throw "No frequency rate estimate available"; }
  virtual CMACK_MSGTYPE frequencyOffsetEstimate() { throw "No frequency offset estimate available"; }
  virtual CMACK_MSGTYPE objectiveFunctionValue() { return Qhat; }

  const int L; //total number of pilots
  const int M; //M for M-PSK
  const unsigned int absD; //number of data symbols
  const unsigned int absP; //number of data symbols
    
protected:
  
  inline int sigma(int k) { return z[k % absD].i; }
  
  //pointers
  const std::vector<int>* Dptr; //data positions
  const std::vector<int>* Pptr; //pilot positions
  const std::vector<mycomplex>* pptr; //array of pilot symbols  

  //estimator output variables
  mycomplex chat;

  //working memory
  std::vector<IndexedReal> z;
  std::vector<std::complex<CMACK_MSGTYPE> > g;

  //useful variables
  const CMACK_MSGTYPE w; // = 2.0*scala.math.Pi/M
  const mycomplex eta; // = exp(jw)
  mycomplex nu; // = eta - 1

 private:
  CMACK_MSGTYPE Qhat; //value of objective function
  
};



class CoherentMackenthunWithDoppler : public CoherentMackenthun {

public:
  CoherentMackenthunWithDoppler(
		     const std::vector<int>* D,
		     const std::vector<int>* P,
		     const std::vector<mycomplex>* p,
		     const int M,
		     const CMACK_MSGTYPE fmin,
		     const CMACK_MSGTYPE fmax,
		     const CMACK_MSGTYPE T); 
  
  virtual void estimate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p);
  virtual inline CMACK_MSGTYPE frequencyOffsetEstimate() { return fhat; }
  virtual inline CMACK_MSGTYPE objectiveFunctionValue() { return Qhat; }

 protected:
  CMACK_MSGTYPE       fhat; //frequency estimate
  const CMACK_MSGTYPE Ts; //symbol period
  CMACK_MSGTYPE fstep;
  const CMACK_MSGTYPE fmin;
  const CMACK_MSGTYPE fmax;
  std::vector<std::complex<CMACK_MSGTYPE> > yf; //working memory

  ///Refine the Doppler estimate using Newton Raphson
  virtual void refine(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p);
  virtual void hardDecisionsAndDerotate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p);
  std::vector<mycomplex> X; //working memory for Newton Raphson

private:
  CMACK_MSGTYPE Qhat; //value of objective function
  
};


class CoherentMackenthunWithDopplerAndDopplerRate : public CoherentMackenthunWithDoppler {

 public:
  CoherentMackenthunWithDopplerAndDopplerRate(
		     const std::vector<int>* D,
		     const std::vector<int>* P,
		     const std::vector<mycomplex>* p,
		     const int M,
		     const CMACK_MSGTYPE fmin,
		     const CMACK_MSGTYPE fmax,
		     const CMACK_MSGTYPE frmin,
		     const CMACK_MSGTYPE frmax,
		     const CMACK_MSGTYPE T); 
  
  virtual void estimate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p);
  virtual inline CMACK_MSGTYPE frequencyRateEstimate() { return frhat; }
  virtual inline CMACK_MSGTYPE objectiveFunctionValue() { return Qhat; }
  virtual void refine(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p);
  virtual void hardDecisionsAndDerotate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p);

 protected:
  CMACK_MSGTYPE       frhat; //frequency rate estimate
  CMACK_MSGTYPE frstep;
  const CMACK_MSGTYPE frmin;
  const CMACK_MSGTYPE frmax;
  ///Refine the Doppler estimate using Newton Raphson
  //virtual void refine(const std::vector<complex>& y, const std::vector<complex>& p);

private:
  CMACK_MSGTYPE Qhat; //value of objective function
};

#endif
