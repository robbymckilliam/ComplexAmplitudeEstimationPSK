//
//  FrequencyEstimator.cpp
//  FrequencyEstimator
//
//  Created by Andre Pollok on 07/02/2013.
//  Copyright (c) 2013 ITR, University of South Australia. All rights reserved.
//

#include <iostream>
# include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <limits>
#include <complex>
#include <vector>
#include <algorithm>
using namespace std;
#include <fftw3.h>

#ifdef __ARM_NEON__
  #include "arm_neon.h"
  #include "neon_mathfun.h"
  #define applyPhasePolynomial applyPhasePolynomialNEON
  // FFTW only available in single-precision (float) on ARM
  #warning "Compiling with float version of FFTW..."
  #define fftw_destroy_plan fftwf_destroy_plan
  #define fftw_plan fftwf_plan
  #define fftw_plan_dft_1d fftwf_plan_dft_1d
  #define fftw_complex fftwf_complex
  #define fftw_cleanup fftwf_cleanup
  #define fftw_execute fftwf_execute
#else
  #define applyPhasePolynomial applyPhasePolynomialC
#endif

#ifndef FFTWFLOAT
  #define FFTWFLOAT 0
#endif
#if FFTWFLOAT==1
  #warning "Compiling with float version of FFTW..."
  #define fftw_destroy_plan fftwf_destroy_plan
  #define fftw_plan fftwf_plan
  #define fftw_plan_dft_1d fftwf_plan_dft_1d
  #define fftw_complex fftwf_complex
  #define fftw_cleanup fftwf_cleanup
  #define fftw_execute fftwf_execute
#endif

#include "FrequencyEstimator.h"

#ifndef M_PI
#define M_PI (3.141592653589793)
#endif 
#ifndef M_PI_4
#define M_PI_4 (M_PI/4.0)
#endif 

//%%% FrequencyEstimator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/// \brief  Standard constructor.
///
FrequencyEstimator::FrequencyEstimator()
{
}

//-----------------------------------------------------------------------------

/// \brief  Standard destructor.
///
FrequencyEstimator::~FrequencyEstimator(void)
{  
}

//=============================================================================

/// \brief  Standard constructor.
/// \param  Nfft    FFT size.
/// \param  Ts      Sample period.
/// \param  tauref  Time reference point for sample times n*Ts-tauref.
///                 Default value is tauref = 0.
///
FrequencyEstimatorLeastSquares::FrequencyEstimatorLeastSquares(const unsigned int Nfft,
                                                               const VALTYPE Ts,
                                                               const VALTYPE tauref) :
  FrequencyEstimator(),
  m_Nfft(Nfft),
  m_Ts(Ts),
  m_tauref(tauref)
{
  // allocate memory for FFT and initialise with zeros
  m_yvec.resize(m_Nfft, complexT(0.0,0.0));
  m_Yvec.resize(m_Nfft, complexT(0.0,0.0));
  m_I.resize(m_Nfft);
  m_X.resize(m_Nfft);
  
  // set up vector with sample times
  m_tSamp.resize(m_Nfft);
  for (unsigned int i=0; i<m_Nfft; i++) {
    m_tSamp[i] = i*m_Ts-m_tauref;
  }

  // create FFTW plan
  m_fftwplan = fftw_plan_dft_1d(m_Nfft, reinterpret_cast<fftw_complex*>(&(m_yvec[0])), reinterpret_cast<fftw_complex*>(&(m_Yvec[0])), FFTW_FORWARD, FFTW_ESTIMATE);
  
  // initialise estimates
  m_fhat = 0;
  m_cgainhat = complexT(0,0);
  m_funVal = -1;
}

//-----------------------------------------------------------------------------

/// \brief  Standard destructor.
///
FrequencyEstimatorLeastSquares::~FrequencyEstimatorLeastSquares(void)
{
  fftw_destroy_plan(m_fftwplan);
  fftw_cleanup();
}

//-----------------------------------------------------------------------------

/// \brief  Runs the frequency offset estimator.
/// \param  rxSamp      vector of received samples
/// \param  txSamp      vector of hypothesised transmit samples
///
void FrequencyEstimatorLeastSquares::estimate(const std::vector<complexT>& rxSamp, const std::vector<complexT>& txSamp)
{
  // number of received samples
  m_numSamp = rxSamp.size();
  
  checkVectorSize(rxSamp, txSamp);
  
  // "demodulate" rx samples using the hypothesised tx samples
  // and compute the normalisation factor for FFT and refinement
  VALTYPE normfactorFFT;
  if (txSamp.size() == 1) { // implicitly assume that txSamp is the all-one vector (i.e. tone at DC)
    for (unsigned int i=0; i<m_numSamp; i++) {
      m_yvec[i] = rxSamp[i];
    }
    // compute normalisation of least-squares estimator (used in refinement)
    normfactorFFT = std::min(m_numSamp,m_Nfft);
    m_normfactor = m_numSamp;
  } else {
    normfactorFFT = 0.0;
    for (unsigned int i=0; i<std::min(m_numSamp,m_Nfft); i++) {
      m_yvec[i] = rxSamp[i] * std::conj(txSamp[i]);
      // compute normalisation of least-squares estimator
      normfactorFFT += std::norm(txSamp[i]);
    }
    m_normfactor = normfactorFFT;
    for (unsigned int i=std::min(m_numSamp,m_Nfft); i<m_numSamp; i++) {
      m_yvec[i] = rxSamp[i] * std::conj(txSamp[i]);
      // compute normalisation of least-squares estimator
      m_normfactor += std::norm(txSamp[i]);
    }
  }
  normfactorFFT = 1.0/normfactorFFT;
  m_normfactor  = 1.0/m_normfactor;
  
  // compute FFT of m_yvec using FFTW3 library:
  fftw_execute(m_fftwplan);

  // compute squared magnitude of FFT of m_yvec
  for (unsigned int i=0; i<m_Nfft; i++) {
    m_I[i] = std::norm(m_Yvec[i]);
  }
  
  // find maximum of objective function and corresponding FFT bin
  std::vector<VALTYPE>::iterator itr_max = std::max_element(m_I.begin(), m_I.end());
  unsigned int idxmax = std::distance(m_I.begin(), itr_max);
  // maximum objective function value
  m_funVal = *itr_max*normfactorFFT*normfactorFFT;
  // complex gain estimate
  m_cgainhat = m_Yvec[idxmax]*normfactorFFT;
  //cout << idxmax << ": " << *itr_max << "\n";
  
  // compute frequency estimate
  m_fhat = (VALTYPE)idxmax/(VALTYPE)m_Nfft/m_Ts;
  // handle negative frequency offsets
  if (m_fhat*2*m_Ts > 1) {
    m_fhat = m_fhat - 1/m_Ts;
  }
  
//  cout << endl << "before refinement:" << endl << "\t" <<
//    "cgainhat=" << m_cgainhat << endl << "\t" <<
//    "fhat=" << m_fhat << endl << "\t" <<
//    "funVal=" << m_funVal << endl;
}

//-----------------------------------------------------------------------------

/// \brief  Runs the frequency offset estimator.
/// \param  rxSamp      vector of received samples
///
/// The transmit signal is assumed to be an all-one vector (i.e. tone at DC).
///
void FrequencyEstimatorLeastSquares::estimate(const std::vector<complexT>& rxSamp)
{
  // no vector of hypothesised transmit samples given -> assume all-one vector (i.e. tone at DC)
  std::vector<complexT> txSamp;
  txSamp.push_back(complexT(1.0,0.0));
  
  // call overloaded estimate()
  estimate(rxSamp, txSamp);
}

//-----------------------------------------------------------------------------

/// \brief  Checks for the correct size of the input vectors.
/// \param  rxSamp      vector of received samples
/// \param  txSamp      vector of hypothesised transmit samples
///
/// If the input vectors are longer than the FFT size, some vector-valued member variables are resized.
/// Note that this requires a reallocation of memory.
/// If the input vectors are shorter than the FFT size, the function ensures that the data is zero padded
/// appropriately before carrying out the FFT.
///
void FrequencyEstimatorLeastSquares::checkVectorSize(const std::vector<complexT>& rxSamp, const std::vector<complexT>& txSamp)
{
  if ( (txSamp.size() != 1) & (txSamp.size() != rxSamp.size()) ) {
    std::ostringstream str;
    str << "Exception in function FrequencyEstimator::estimate. The vector txSamp must have a length of either 1 or equal to that of rxSamp";
    throw str.str();
  }
  
  // check input vector lengths compared to FFT size
  if (m_numSamp > m_yvec.size()) { // input vectors longer than FFT size -> need to resize
    unsigned int numSampOLD = m_yvec.size();
    m_tSamp.resize(m_numSamp);
    m_X.resize(m_numSamp);
    for (unsigned int i=numSampOLD; i<m_numSamp; i++) {
      m_tSamp[i] = i*m_Ts-m_tauref;
    }
    m_yvec.resize(m_numSamp);
    cout << "Input vectors are longer than specified FFT size. Need to re-create FFTW plan." << endl;
    m_fftwplan = fftw_plan_dft_1d(m_Nfft, reinterpret_cast<fftw_complex*>(&(m_yvec[0])), reinterpret_cast<fftw_complex*>(&(m_Yvec[0])), FFTW_FORWARD, FFTW_ESTIMATE);
  } else if (m_numSamp < m_Nfft) { // input vectors shorter than FFT size -> zero pad
    for (unsigned int i=m_numSamp; i<m_Nfft; i++) {
      m_yvec[i] = complexT(0.0,0.0);
    }
  }
}

//-----------------------------------------------------------------------------

/// \brief  Refine estimates using the Newton Raphson optimiser.
/// \param  tol       Tolerance for optimiser.
///                   Default value is tol=1e-10.
/// \param  max_iter  Maximum number of iterations.
///                   Default value is max_iter=10.
/// \param  alpha     alpha=0: max log(periodogram), alpha=1: max standard periodogram.
///                   Default value is alpha=0.
///
/// Refine the frequency and complex gain estimate using the Newton Raphson optimiser.
/// For details, see Andre Pollok's ASRP workbook, page 34-35.
///
void FrequencyEstimatorLeastSquares::refineEstimate(const VALTYPE tol,
                                                    const unsigned int max_iter,
                                                    const VALTYPE alpha)
{
  complexT  Y, dYdf, dY2df2; // Y and its 1st and 2nd derivative
  VALTYPE   I, dIdf, dI2df2; // I and its 1st and 2nd derivative
  VALTYPE   G, H; // Hessian
  
  const VALTYPE normfactor2pi   = 2.0*M_PI*m_normfactor;
  const VALTYPE normfactor4pi2  = 4.0*M_PI*M_PI*m_normfactor;

  // initialise Newton's method
  VALTYPE fhat = m_fhat; // initialise frequency estimate
  // handle negative frequency offsets
  if (fhat < 0)
    fhat = fhat + 1/m_Ts;
  VALTYPE newtonstep = tol + 1.0; // something bigger than the tolerance to start with
  unsigned int numiters = 0; // number of iterations performed
  
  // run Newton's method
  while( (std::fabs(newtonstep) > tol) && (numiters < max_iter) ){
    Y       = complexT(0,0);
    dYdf    = complexT(0,0);
    dY2df2  = complexT(0,0);
    applyPhasePolynomial(m_yvec, m_X, m_numSamp, 0, -fhat, 0);
    for(unsigned int i = 0; i < m_numSamp; i++) {
      Y       +=                          m_X[i]; // compute Y
      dYdf    += m_tSamp[i]             * m_X[i]; // compute 1st derivative of Y
      dY2df2  += m_tSamp[i]*m_tSamp[i]  * m_X[i]; // compute 2nd derivative of Y
    }
    // normalise
    Y       = complexT(m_normfactor, 0.0)     * Y;
    dYdf    = complexT(0.0, -normfactor2pi)   * dYdf;
    dY2df2  = complexT(-normfactor4pi2, 0.0)  * dY2df2;
    
    I       = std::norm(Y); // compute I
    dIdf    = 2.0*std::real( dYdf   * std::conj(Y) ); // 1st derivative of I
    dI2df2  = 2.0*std::real( dY2df2 * std::conj(Y) ) + 2.0*std::norm(dYdf); // 2nd derivative of I
    G = dIdf*dIdf;
    H = dI2df2;
    
    // Newton iterate with monotonic function kappa(I) = I^alpha.
    // For alpha=1, this is the standard Newton iterate inv(H)*[dIdf; dIdfr].
    newtonstep = dIdf / ((alpha-1)/I*G + H);
    fhat = fhat - newtonstep; // update fhat
    
    numiters += 1;
  }
  
  // final frequency estimate
  m_fhat = fhat;
  // handle negative frequency offsets
  if (m_fhat*2*m_Ts > 1)
    m_fhat = m_fhat - 1/m_Ts;

  // final complex gain estimate
  m_cgainhat = Y;
  
  // final objective function value
  m_funVal = std::norm(Y);
  
//  cout << endl << "after refinement:" << endl << "\t" <<
//    "cgainhat=" << m_cgainhat << endl << "\t" <<
//    "fhat=" << m_fhat << endl << "\t" <<
//    "funVal=" << m_funVal << endl;
}

//-----------------------------------------------------------------------------

/// \brief  Multiplies complex input vector with a complex exponential.
/// \param  in        Complex input vector.
/// \param  out       Complex output vector.
/// \param  numSamp   Number of samples to be processed (can be smaller than the input vector size).
/// \param  phi       Zeroth order coefficient (phase) of polynomial phase signal.
/// \param  f         First order coefficient (frequency) of polynomial phase signal.
/// \param  fr        Second order coefficient (frequency rate) of polynomial phase signal.
///
/// Multiplies a complex input vector with a complex exponential, whose phase is determined by the
/// specified polynomial phase coefficients. The polynomial phase signal is:
///
///         phi + 2*pi * (f*t + (fr/2)*t^2)
///
void FrequencyEstimatorLeastSquares::applyPhasePolynomialC(const std::vector<complexT>& in,
                                                          std::vector<complexT>& out,
                                                          const unsigned int numSamp,
                                                          const VALTYPE phi,
                                                          const VALTYPE f,
                                                          const VALTYPE fr)
{
	VALTYPE theta;
  const VALTYPE f2pi  = 2.0*M_PI*f;
  const VALTYPE fr2pi = 2.0*M_PI*(0.5*fr);
  
	for(unsigned int i=0; i<numSamp; i++)
	{
		theta = phi + f2pi*m_tSamp[i] + fr2pi*m_tSamp[i]*m_tSamp[i];
    
    // standard C:
    out[i]  = in[i] * std::polar<VALTYPE>(1.0, theta);

    // fast cos/sin from phy_utilities:
		//out[i]  = in[i] * complexT(cos_32_large(theta), sin_32_large(theta));
	}
}

//-----------------------------------------------------------------------------

/// \brief  Multiplies complex input vector with a complex exponential.
/// \param  in        Complex input vector.
/// \param  out       Complex output vector.
/// \param  phi       Zeroth order coefficient (phase) of polynomial phase signal.
/// \param  f         First order coefficient (frequency) of polynomial phase signal.
/// \param  fr        Second order coefficient (frequency rate) of polynomial phase signal.
///
/// Multiplies a complex input vector with a complex exponential, whose phase is determined by the
/// specified polynomial phase coefficients. The polynomial phase signal is:
///
///         phi + 2*pi * (f*t + (fr/2)*t^2)
///
void FrequencyEstimatorLeastSquares::applyPhasePolynomialC(const std::vector<complexT>& in,
                                                          std::vector<complexT>& out,
                                                          const VALTYPE phi,
                                                          const VALTYPE f,
                                                          const VALTYPE fr)
{
  applyPhasePolynomial(in, out, in.size(), phi, f, fr);
}

//-----------------------------------------------------------------------------

#ifdef __ARM_NEON__
/// \brief  Multiplies complex input vector with a complex exponential (ARM NEON implementation).
/// \param  in        Complex input vector.
/// \param  out       Complex output vector.
/// \param  numSamp   Number of samples to be processed (can be smaller than the input vector size).
/// \param  phi       Zeroth order coefficient (phase) of polynomial phase signal.
/// \param  f         First order coefficient (frequency) of polynomial phase signal.
/// \param  fr        Second order coefficient (frequency rate) of polynomial phase signal.
///
/// CAUTION: implementation uses ARM NEON instructions!
///
/// Multiplies a complex input vector with a complex exponential, whose phase is determined by the
/// specified polynomial phase coefficients. The polynomial phase signal is:
///
///         phi + 2*pi * (f*t + (fr/2)*t^2)
///
void FrequencyEstimatorLeastSquares::applyPhasePolynomialNEON(const std::vector<complexT>& in,
                                                              std::vector<complexT>& out,
                                                              const unsigned int numSamp,
                                                              const VALTYPE phi,
                                                              const VALTYPE f,
                                                              const VALTYPE fr)
{
	unsigned int	n;
	
	float32x4x2_t	vin, vout;
	float32x4_t		vt, vt2, vtheta, vsin, vcos;
		
	// precompute a few scalar products
	const VALTYPE f2pi  	= 2*M_PI*f;
	const VALTYPE fr2pi 	= 2*M_PI*(0.5*fr);
	
	for(n=0; n<numSamp; n+=4)
	{
		// load sample times
    vt   = vld1q_f32((float32_t*)&(m_tSamp[n]));
    // squared sample times
		vt2  = vmulq_f32(vt, vt);
		
		// compute theta
		vtheta = vld1q_dup_f32((float32_t*)&phi);
		vtheta = vmlaq_n_f32(vtheta, vt, f2pi);
		vtheta = vmlaq_n_f32(vtheta, vt2, fr2pi);
    
		// compute sin and cos
		sincos_ps(vtheta, &vsin, &vcos);
		
		// load 4 complex values
		vin  = vld2q_f32((float32_t*)&(in[n]));
		
		// compute complex product (val[0] is the real part and val[1] is the imaginary part)
		vout.val[0] = vmulq_f32(vin.val[0], vcos);
		vout.val[0] = vmlsq_f32(vout.val[0], vin.val[1], vsin);
		vout.val[1] = vmulq_f32(vin.val[0], vsin);
		vout.val[1] = vmlaq_f32(vout.val[1], vin.val[1], vcos);
		
		// save 4 complex values
		vst2q_f32((float32_t*)&(out[n]), vout);
  }
}
#endif

//=============================================================================

/// \brief  Standard constructor.
/// \param  Nfft                FFT size.
/// \param  Ts                  Sample period.
/// \param  tauref              Time reference point for sample times n*Ts-tauref.
/// \param  frmin               Lower limit of the frequency rate search range.
/// \param  frmax               Upper limit of the frequency rate search range.
/// \param  search_oversample   Oversampling factor for the frequency rate search.
///                             Default value is tauref = 1.0.
///
FrequencyEstimatorLeastSquaresWithRate::FrequencyEstimatorLeastSquaresWithRate(const unsigned int Nfft,
                                                                               const VALTYPE Ts,
                                                                               const VALTYPE tauref,
                                                                               const VALTYPE frmin,
                                                                               const VALTYPE frmax,
                                                                               const VALTYPE search_oversample) :
  FrequencyEstimatorLeastSquares(Nfft,Ts,tauref),
  m_frmin(frmin),
  m_frmax(frmax),
  m_search_oversample(search_oversample)
{  
  m_frhat = 0;
}

//-----------------------------------------------------------------------------

/// \brief  Standard destructor.
///
FrequencyEstimatorLeastSquaresWithRate::~FrequencyEstimatorLeastSquaresWithRate(void)
{
}

//-----------------------------------------------------------------------------

/// \brief  Runs the frequency and frequency rate estimator.
/// \param  rxSamp      vector of received samples
/// \param  txSamp      vector of hypothesised transmit samples
///
void FrequencyEstimatorLeastSquaresWithRate::estimate(const std::vector<complexT>& rxSamp, const std::vector<complexT>& txSamp)
{
  // number of received samples
  m_numSamp = rxSamp.size();
    
  checkVectorSize(rxSamp, txSamp);
  
  // "demodulate" rx samples using the hypothesised tx samples
  // and compute the normalisation factor for FFT and refinement
  VALTYPE normfactorFFT;
  if (txSamp.size() == 1) { // implicitly assume that txSamp is the all-one vector (i.e. tone at DC)
    for (unsigned int i=0; i<m_numSamp; i++) {
      m_yvec[i] = rxSamp[i];
    }
    // compute normalisation of least-squares estimator (used in refinement)
    normfactorFFT = std::min(m_numSamp,m_Nfft);
    m_normfactor = m_numSamp;
  } else {
    normfactorFFT = 0.0;
    for (unsigned int i=0; i<std::min(m_numSamp,m_Nfft); i++) {
      m_yvec[i] = rxSamp[i] * std::conj(txSamp[i]);
      // compute normalisation of least-squares estimator
      normfactorFFT += std::norm(txSamp[i]);
    }
    m_normfactor = normfactorFFT;
    for (unsigned int i=std::min(m_numSamp,m_Nfft); i<m_numSamp; i++) {
      m_yvec[i] = rxSamp[i] * std::conj(txSamp[i]);
      // compute normalisation of least-squares estimator
      m_normfactor += std::norm(txSamp[i]);
    }
  }
  normfactorFFT = 1.0/normfactorFFT;
  m_normfactor  = 1.0/m_normfactor;
    
  setFreqRateSearchGrid();
  VALTYPE fr = m_frmin;
  std::vector<complexT> yvec(m_yvec); // create backup of m_yvec
  for (unsigned int i=0; i<m_frtrials; i++) {
    //cout << "trialing fr=" << fr << endl;
    
    // compensate by trial frequency rate of change
    applyPhasePolynomial(yvec, m_yvec, std::min(m_Nfft,m_numSamp), 0, 0, -fr);
    
    // compute FFT of m_yvec using FFTW3 library
    fftw_execute(m_fftwplan);

    // compute squared magnitude of FFT of m_yvec
    for (unsigned int i=0; i<m_Nfft; i++) {
      m_I[i] = std::norm(m_Yvec[i]);
    }
    
    // find maximum of objective function and corresponding FFT bin
    std::vector<VALTYPE>::iterator itr_max = std::max_element(m_I.begin(), m_I.end());
    unsigned int idxmax = std::distance(m_I.begin(), itr_max);
    //cout << "fr=" << fr << ": " << idxmax << ", " << *itr_max << endl;
    
    if (m_funVal < *itr_max) {
      // maximum objective function value
      m_funVal = *itr_max;
      // get corresponding complex gain, frequency and frequency rate estimates
      m_cgainhat = m_Yvec[idxmax];
      m_fhat  = (VALTYPE)idxmax/(VALTYPE)m_Nfft/m_Ts;
      m_frhat = fr;
    }
    
    fr += m_frstep;
  }
  
  // normalisation
  m_funVal    *= normfactorFFT*normfactorFFT;
  m_cgainhat  *= normfactorFFT;
    
  // restore original m_yvec
  m_yvec.swap(yvec);
    
  // handle negative frequency offsets
  if (m_fhat*2*m_Ts > 1) {
    m_fhat = m_fhat - 1/m_Ts;
  }

//  cout << endl << "before refinement:" << endl << "\t" <<
//    "cgainhat=" << m_cgainhat << endl << "\t" <<
//    "fhat=" << m_fhat << endl << "\t" <<
//    "frhat=" << m_frhat << endl << "\t" <<
//    "funVal=" << m_funVal << endl;
}

//-----------------------------------------------------------------------------

/// \brief  Sets grid point spacing and number of trials for frequency rate search.
///
void FrequencyEstimatorLeastSquaresWithRate::setFreqRateSearchGrid(void)
{
  m_frtrials = (unsigned int)ceil(std::abs(m_frmax - m_frmin) * m_Ts*m_Ts*m_numSamp*m_numSamp*m_search_oversample) + 1;
  
  if (m_frtrials == 0)
    m_frtrials = 1; // at least one trial
  
  m_frstep = std::abs(m_frmax - m_frmin)/(m_frtrials-1);
  
  if (m_frstep == 0.0)
    m_frstep = 1.0; // prevent potential infinite loop
  
  //std::cout << "number of trials = " << frtrials << ", stepsize = " << m_frstep << ", last trial fr = " << m_frmin + (frtrials-1)*m_frstep << std::endl;
}

//-----------------------------------------------------------------------------

/// \brief  Refine estimates using the Newton Raphson optimiser.
/// \param  tol       Tolerance for optimiser.
///                   Default value is tol=1e-10.
/// \param  max_iter  Maximum number of iterations.
///                   Default value is max_iter=10.
/// \param  alpha     alpha=0: max log(periodogram), alpha=1: max standard periodogram.
///                   Default value is alpha=0.
///
/// Refine the frequency and complex gain estimate using the Newton Raphson optimiser.
/// For details, see Andre Pollok's ASRP workbook, page 34-35.
///
void FrequencyEstimatorLeastSquaresWithRate::refineEstimate(const VALTYPE tol,
                                                    const unsigned int max_iter,
                                                    const VALTYPE alpha)
{
  // Y and its 1st and 2nd derivative
  complexT  Y,
            dYdf, dYdfr,
            dY2df2, dY2dfr2, dY2dfrdf;
  // I and its 1st and 2nd derivative
  VALTYPE   I, Iinv,
            dIdf, dIdfr,
            dI2df2, dI2dfr2, dI2dfrdf;
  // Hessian matrix and its reciprocal determinant and matrix inverse
  VALTYPE H11, H12, H21, H22, detinv, Hinv11, Hinv12, Hinv21, Hinv22;
  
  const VALTYPE normfactor2pi   = 2.0*M_PI*m_normfactor;
  const VALTYPE normfactor4pi2  = 4.0*M_PI*M_PI*m_normfactor;
  
  // initialise Newton's method
  VALTYPE fhat   = m_fhat; // initialise frequency estimate
  // handle negative frequency offsets
  if (fhat < 0)
    fhat = fhat + 1/m_Ts;
  VALTYPE etahat = 0.5*m_frhat; // initialise frequency rate estimate (factor 2 due to frequency rate = 2nd derivative of (f*t + eta*t^2) = 2*eta)
  VALTYPE newtonstepf, newtonstepeta;
  VALTYPE newtonstep = tol + 1.0; // something bigger than the tolerance to start with
  unsigned int numiters = 0; // number of iterations performed
  // run Newton's method
  while( (std::fabs(newtonstep) > tol) && (numiters < max_iter) ){
    Y         = complexT(0,0);
    dYdf      = complexT(0,0);
    dYdfr     = complexT(0,0);
    dY2df2    = complexT(0,0);
    dY2dfr2   = complexT(0,0);
    dY2dfrdf  = complexT(0,0);
    applyPhasePolynomial(m_yvec, m_X, m_numSamp, 0, -fhat, -2.0*etahat);
    for(unsigned int i = 0; i < m_numSamp; i++) {
      Y         +=                                      m_X[i]; // compute Y
      dYdf      += m_tSamp[i]                         * m_X[i]; // 1st derivative of Y with respect to f
      dYdfr     += m_tSamp[i]*m_tSamp[i]              * m_X[i]; // 1st derivative of Y with respect to fr
      dY2df2    += m_tSamp[i]*m_tSamp[i]              * m_X[i]; // 2nd derivative of Y with respect to f
      dY2dfr2   += std::pow(m_tSamp[i], (VALTYPE)4.0) * m_X[i]; // 2nd derivative of Y with respect to fr
      dY2dfrdf  += std::pow(m_tSamp[i], (VALTYPE)3.0) * m_X[i]; // derivative of Y with respect to fr and f
    }
    // normalise
    Y         = complexT(m_normfactor, 0.0)     * Y;
    dYdf      = complexT(0.0, -normfactor2pi)   * dYdf;
    dYdfr     = complexT(0.0, -normfactor2pi)   * dYdfr;
    dY2df2    = complexT(-normfactor4pi2, 0.0)  * dY2df2;
    dY2dfr2   = complexT(-normfactor4pi2, 0.0)  * dY2dfr2;
    dY2dfrdf  = complexT(-normfactor4pi2, 0.0)  * dY2dfrdf;
    
    I         = std::norm(Y); //compute I
    dIdf      = 2.0*std::real( dYdf     * std::conj(Y) ); // 1st derivative of I with respect to f
    dI2df2    = 2.0*std::real( dY2df2   * std::conj(Y) )  + 2.0*std::norm(dYdf); // 2nd derivative with respect to f
    dIdfr     = 2.0*std::real( dYdfr    * std::conj(Y) ); // 1st derivative with respect to fr
    dI2dfr2   = 2.0*std::real( dY2dfr2  * std::conj(Y) )  + 2.0*std::norm(dYdfr); // 2nd derivative with respect to fr
    dI2dfrdf  = 2.0*std::real( dY2dfrdf * std::conj(Y)    + dYdf * std::conj(dYdfr) ); // derivative with respect to fr and f
    
    // compute 2x2 Hessian matrix and its inverse
    Iinv = (alpha-1.0)/I;
    H11 = Iinv * dIdf*dIdf   + dI2df2;
    H12 = Iinv * dIdf*dIdfr  + dI2dfrdf;
    H21 = H12;
    H22 = Iinv * dIdfr*dIdfr + dI2dfr2;
    detinv = 1.0/(H11*H22 - H12*H21); // determinant
    Hinv11 =  H22*detinv;
    Hinv12 = -H12*detinv;
    Hinv21 = -H21*detinv;
    Hinv22 =  H11*detinv;
    
    // Newton iterate with monotonic function kappa(I) = I^alpha.
    // For alpha=1, this is the standard Newton iterate inv(H)*[dIdf; dIdfr].
    newtonstepf    = Hinv11*dIdf + Hinv12*dIdfr;
    newtonstepeta  = Hinv21*dIdf + Hinv22*dIdfr;
    fhat   = fhat   - newtonstepf;   // update fhat
    etahat = etahat - newtonstepeta; // update frhat
    
    //step size for termination
    newtonstep = std::fabs(newtonstepf) + std::fabs(newtonstepeta);
 
    numiters += 1;
  }
  
  // final frequency and frequency rate estimates
  m_fhat  = fhat;
  // handle negative frequency offsets
  if (m_fhat*2*m_Ts > 1)
    m_fhat = m_fhat - 1/m_Ts;
  m_frhat = 2*etahat; // factor 2 due to frequency rate = 2nd derivative of (f*t + eta*t^2) = 2*eta

  // final complex gain estimate
  m_cgainhat = Y;
  
  // final objective function value
  m_funVal = std::norm(Y);
  
//  cout << endl << "after refinement:" << endl << "\t" <<
//    "cgainhat=" << m_cgainhat << endl << "\t" <<
//    "fhat=" << m_fhat << endl << "\t" <<
//    "frhat=" << m_frhat << endl << "\t" <<
//    "funVal=" << m_funVal << endl;
}

//-----------------------------------------------------------------------------















