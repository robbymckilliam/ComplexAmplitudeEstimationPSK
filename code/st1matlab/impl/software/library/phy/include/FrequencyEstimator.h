//
//  FrequencyEstimator.h
//  FrequencyEstimator
//
//  Created by Andre Pollok on 07/02/2013.
//  Copyright (c) 2013 ITR, University of South Australia. All rights reserved.
//

#ifndef FrequencyEstimator_FrequencyEstimator_h
#define FrequencyEstimator_FrequencyEstimator_h

#ifndef VALTYPE
#define VALTYPE double
#endif
#ifndef MSGTYPE
#define MSGTYPE VALTYPE
#endif

#include <fftw3.h>

typedef std::complex<VALTYPE> complexT;

class FrequencyEstimator
{
  
public:
  FrequencyEstimator();
  ~FrequencyEstimator(void);
  
  virtual VALTYPE frequencyEstimate() = 0;
  virtual VALTYPE frequencyRateEstimate() = 0;
  virtual complexT complexGainEstimate() = 0;
  virtual VALTYPE objectiveFunctionValue() = 0;
  
  virtual void estimate(const std::vector<complexT>& rxSamp) = 0;
  virtual void estimate(const std::vector<complexT>& rxSamp, const std::vector<complexT>& txSamp) = 0;
  
  virtual void refineEstimate(const VALTYPE tol=1e-10,
                              const unsigned int max_iter=10,
                              const VALTYPE alpha=0.0) = 0;
  
protected:
  
};

class FrequencyEstimatorLeastSquares : public FrequencyEstimator
{
  
public:
  FrequencyEstimatorLeastSquares(const unsigned int Nfft,
                                 const VALTYPE Ts,
                                 const VALTYPE tauref);
  ~FrequencyEstimatorLeastSquares(void);
  
  virtual inline VALTYPE frequencyEstimate() { return m_fhat; };
  virtual inline VALTYPE frequencyRateEstimate() { throw "No frequency rate estimate available"; };
  virtual inline complexT complexGainEstimate() { return m_cgainhat; };
  virtual inline VALTYPE objectiveFunctionValue() { return m_funVal; };
  
  virtual void estimate(const std::vector<complexT>& rxSamp);
  virtual void estimate(const std::vector<complexT>& rxSamp,
                        const std::vector<complexT>& txSamp);
  
  virtual void refineEstimate(const VALTYPE tol=1e-10,
                              const unsigned int max_iter=10,
                              const VALTYPE alpha=0.0);
    
protected:
  
  virtual void checkVectorSize(const std::vector<complexT>& rxSamp,
                               const std::vector<complexT>& txSamp);

  virtual void applyPhasePolynomialC(const std::vector<complexT>& in,
                                     std::vector<complexT>& out,
                                     const unsigned int numSamp,
                                     const VALTYPE phi,
                                     const VALTYPE f,
                                     const VALTYPE fr);
  virtual void applyPhasePolynomialC(const std::vector<complexT>& in,
                                     std::vector<complexT>& out,
                                     const VALTYPE phi,
                                     const VALTYPE f,
                                     const VALTYPE fr);
#ifdef __ARM_NEON__
  virtual void applyPhasePolynomialNEON(const std::vector<complexT>& in,
                                        std::vector<complexT>& out,
                                        const unsigned int numSamp,
                                        const VALTYPE phi,
                                        const VALTYPE f,
                                        const VALTYPE fr);
#endif

  // estimates of complex gain and frequency
  complexT m_cgainhat;
  VALTYPE m_fhat;
  // final objective function value
  VALTYPE m_funVal;
  
  // FFT size and number of samples
  const unsigned int m_Nfft;
  unsigned int m_numSamp;
  
  // sample period and time reference point
  const VALTYPE m_Ts;
  const VALTYPE m_tauref;
  std::vector<VALTYPE> m_tSamp;

  // normalisation of least-squares estimator
  VALTYPE m_normfactor;
  
  // input, output and plan for FFTW
  std::vector<complexT> m_yvec;
  std::vector<complexT> m_Yvec;
  // squared absolute value of FFT output
  std::vector<VALTYPE> m_I;
  std::vector<complexT> m_X;
  fftw_plan m_fftwplan;
};

class FrequencyEstimatorLeastSquaresWithRate : public FrequencyEstimatorLeastSquares
{
  
public:
  FrequencyEstimatorLeastSquaresWithRate(const unsigned int Nfft,
                                         const VALTYPE Ts,
                                         const VALTYPE tauref,
                                         const VALTYPE frmin,
                                         const VALTYPE frmax,
                                         const VALTYPE search_oversample=1.0);
  ~FrequencyEstimatorLeastSquaresWithRate(void);
  
  virtual inline VALTYPE frequencyRateEstimate() { return m_frhat; };
    
  virtual void estimate(const std::vector<complexT>& rxSamp, const std::vector<complexT>& txSamp);
  
  virtual void refineEstimate(const VALTYPE tol=1e-10,
                              const unsigned int max_iter=10,
                              const VALTYPE alpha=0.0);
  
protected:
  
  virtual void setFreqRateSearchGrid(void);
 
  // estimate of frequency rate of change
  VALTYPE m_frhat;
    
  // bounds, oversampling factor and step size for frequency rate search range
  const VALTYPE m_frmin;
  const VALTYPE m_frmax;
  const VALTYPE m_search_oversample;
  unsigned int  m_frtrials;
  VALTYPE       m_frstep;
};

#endif

