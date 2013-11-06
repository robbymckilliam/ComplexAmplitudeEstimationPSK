//
//  TimingEstimator.cpp
//  TimingEstimator
//
//  Created by Andre Pollok on 31/05/2012.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#include <iostream>
# include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <limits>
#include <complex>
using namespace std;
#include "TimingEstimator.h"


#ifndef M_PI
#define M_PI (3.141592653589793)
#endif 
#ifndef M_PI_4
#define M_PI_4 (M_PI/4.0)
#endif 

/// \brief	returns minimum of a and b
///
#define MIN(a,b)			(((a) < (b)) ? (a) : (b))

/// \brief	returns maximum of a and b
///
#define MAX(a,b)			(((a) > (b)) ? (a) : (b))

//%%% PulseShape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/// \brief	Standard constructor.
/// \param  Tsymb       symbol period
/// \param  OS          oversampling factor
/// \param  duration    total duration of the pulse
///
PulseShape::PulseShape(VALTYPE Tsymb, VALTYPE OS, VALTYPE duration) 
  : m_Tsymb(Tsymb),
    m_OS(OS),
    m_duration(duration)
{
}

//-----------------------------------------------------------------------------

/// \brief	Returns the symbol period.
///
VALTYPE PulseShape::getTsymb()
{
  return m_Tsymb;
}

//-----------------------------------------------------------------------------

/// \brief	Returns the oversampling factor.
///
VALTYPE PulseShape::getOS()
{
  return m_OS;
}

//-----------------------------------------------------------------------------

/// \brief	Returns the duration of the pulse.
///
VALTYPE PulseShape::getduration()
{
  return m_duration;
}

//=============================================================================

/// \brief	Standard constructor for root raised cosine pulse object.
/// \param  rolloff     roll off factor between 0 and 1
/// \param  Tsymb       symbol period
/// \param  OS          oversampling factor
/// \param  duration    total duration of the pulse
///
RootRaisedCosine::RootRaisedCosine(VALTYPE rolloff, VALTYPE Tsymb, VALTYPE OS, VALTYPE duration)
  : PulseShape(Tsymb, OS, duration),
    m_rolloff(rolloff)
{
  m_fc = 1/(2*m_Tsymb);       // cut off frequency
  m_fs = m_OS/m_Tsymb;    // sampling frequency
  m_tol = sqrt(std::numeric_limits<VALTYPE>::epsilon());
}

//-----------------------------------------------------------------------------

/// \brief	Returns the value of the root raised cosine pulse at the given time.
///         The implementation based on MATLAB's firrcos.
/// \param  time        time at which the pulse is to be evaluated
///
VALTYPE RootRaisedCosine::pulseshape(VALTYPE time)
{  
  VALTYPE rrc;
    
  if (abs(time) > m_duration) {
    rrc = 0;
  }
  else if (abs(time) == 0) {
    rrc = -sqrt(2*m_fc) / (M_PI*m_fs) * (M_PI*(m_rolloff-1) - 4*m_rolloff);
  }
  else if (abs(abs(8*m_rolloff*m_fc*time) - 1.0) < m_tol) {
    rrc = sqrt(2*m_fc) / (2*M_PI*m_fs)
    * (M_PI*(m_rolloff+1)  * sin(M_PI_4*(m_rolloff+1)/m_rolloff)
       - 4*m_rolloff * sin(M_PI_4*(m_rolloff-1)/m_rolloff)
       + M_PI*(m_rolloff-1) * cos(M_PI_4*(m_rolloff-1)/m_rolloff)
       );
  }
  else {
    rrc = -4*m_rolloff/m_fs
    * ( cos((1+m_rolloff)*2*M_PI*m_fc*time) + sin((1-m_rolloff)*2*M_PI*m_fc*time) / (8*m_rolloff*m_fc*time) )
    / (M_PI * sqrt(1/(2*m_fc)) * (pow(8*m_rolloff*m_fc*time,2) - 1));
  }
  
  // scale pulse
  rrc *= sqrt(2*m_fc);
  rrc /= sqrt(m_Tsymb/m_OS/m_OS);
  
  return rrc;
}

//%%% TimingEstimator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/// \brief  Standard constructor.
///
TimingEstimator::TimingEstimator(VALTYPE          taumin,
                                 VALTYPE          taumax,
                                 PulseShape*      pulseshape,
                                 unsigned int*    D,
                                 unsigned int     numData,
                                 complexT*        pilot,
                                 unsigned int*    P,
                                 unsigned int     numPilots,
                                 unsigned int     a)
{
  unsigned int k, i;
  
  // initialisation
  m_pulseshape = pulseshape;
  m_T = m_pulseshape->getTsymb();
  m_q = m_pulseshape->getOS();
  m_p = 1;
  m_Ts = m_T/m_q;
  m_a = a;
  m_c = m_a*m_p;
  m_Delta = m_Ts/m_c;
  
  m_K = floor((taumax-taumin)/m_Delta)+1;
  m_taugrid = new VALTYPE [m_K];
  m_taugrid[0] = taumin;
  for (k=1; k<m_K; k++) {
    m_taugrid[k] = m_taugrid[k-1] + m_Delta;
  }
  m_gk_tau = taumin - m_Delta; // set to value outside of m_taugrid to ensure m_gk gets computed the first time
  m_bk_tau = taumin - m_Delta;
  m_Yk_tau = taumin - m_Delta;
  m_Zk_tau = taumin - m_Delta;
  
  m_numData   = numData;
  m_numPilots = numPilots;
  m_L = m_numData + m_numPilots;
  
  m_D    = D;
  m_P    = P;
  m_pilot = pilot;
  m_minD = D[0];
  m_maxD = D[numData-1];
  m_minP = P[0];
  m_maxP = P[numPilots-1];
  m_minI = MIN(m_minD,m_minP);
  m_maxI = MAX(m_maxD,m_maxP);
  
  // get a zero filled reflected D indicator
  if (m_maxD-m_minD+1 == m_numData) {
    // contiguous data (reflected D indicator not needed)
    m_hD = NULL;
  } else {
    // non-contiguous data
    m_hD = new bool [m_a*m_q*(m_maxD-m_minD)+1];
    for (i=0; i<m_a*m_q*(m_maxD-m_minD)+1; i++) {
      m_hD[i] = false;
    }
    for (i=0; i<m_numData; i++) {
      m_hD[m_a*m_q*(m_maxD - m_D[i])] = true;
    }
  }
  
  // get a zero filled reflected P indicator
  m_hP = new complexT [m_a*m_q*(m_maxP-m_minP)+1];
  for (i=0; i<m_a*m_q*(m_maxP-m_minP)+1; i++) {
    m_hP[i] = complexT (0,0);
  }
  for (i=0; i<m_numPilots; i++) {
    m_hP[m_a*m_q*(m_maxP - m_P[i])] = conj(m_pilot[i]);
  }
  
  // set default tolerance for the refining step of the estimator (can be changed with setTolerance(tol))
  m_refine_tol = 1e-4;
  
  m_timingEstimate = 0;
  m_funcVal = 0;
    
  m_zk      = NULL;
  m_gk      = NULL;
  m_bk      = NULL;
  m_bkidx   = NULL;
  m_Zk      = NULL;
  m_Yk      = NULL;
}

//-----------------------------------------------------------------------------

/// \brief  Standard destructor.
///
TimingEstimator::~TimingEstimator(void)
{  
  delete[] m_taugrid;

  delete[] m_hD;
  delete[] m_hP;
  
  delete[] m_zk;
  delete[] m_gk;
  delete[] m_bk;
  delete[] m_bkidx;
  delete[] m_Zk;
  delete[] m_Yk;
}

//-----------------------------------------------------------------------------

/// \brief  Runs the timing offset estimator.
/// \param  rxSamp      vector of received samples
/// \param  numSamp     number of samples in the received vector
///
void TimingEstimator::runEstimator(complexT *rxSamp, unsigned int numSamp)
{  
  unsigned int k;
  VALTYPE *SSk = NULL;
  VALTYPE minSSk;
  VALTYPE tauhatcoarse;
    
  m_N = numSamp;
  initzk(rxSamp);

  // compute Yk and Zk
  checkVectorsUpToDate(m_taugrid[0]);
  getYk(m_K, m_taugrid[0]);
  getZk(m_K, m_taugrid[0]);
  // compute sum of squares for all time grid points
  SSk = new VALTYPE [m_K];
  for (k=0; k<m_K; k++) {
    SSk[k] = -pow(m_Zk[k] + abs(m_Yk[k]),2)/m_L;
  }
  
  // find coarse timing offset estimate
  minSSk = SSk[0];
  tauhatcoarse = m_taugrid[0];
  for (k=1; k<m_K; k++) {
    if (SSk[k] < minSSk) {
      minSSk = SSk[k];
      tauhatcoarse = m_taugrid[k];
    }
  }
  
  // hill climb between tauhatcoarse-m_Delta and tauhatcoarse+m_Delta using Brent's method
  refine(tauhatcoarse);
  
  delete[] SSk;
}

//-----------------------------------------------------------------------------

/// \brief  Get zero-filled version of received sequence.
/// \param  rxSamp      vector of received samples
///
void TimingEstimator::initzk(complexT *rxSamp)
{
  unsigned int k;

  if (!m_zk) {
    m_zk = new complexT [m_N*m_c];
  }

  for (k=0; k<m_N*m_c; k++) {
    m_zk[k] = 0;
  }
  for (k=0; k<m_N; k++) {
    m_zk[(k+1)*m_c-1] = rxSamp[k];
  }
}

//-----------------------------------------------------------------------------

/// \brief  Checks if gk, bk, Yk and Zk have been computed for the specified 
///         timing offset. If not, these vectors are deleted.
/// \param  tau       timing offset
///
void TimingEstimator::checkVectorsUpToDate(VALTYPE tau)
{
  initgk(tau);
  if (m_bk_tau != tau) {
    delete[] m_bk;
    m_bk = NULL;
    delete[] m_bkidx;
    m_bkidx = NULL;
    m_bk_tau = tau;
  }
  if (m_Yk_tau != tau) {
    delete[] m_Yk;
    m_Yk = NULL;
    m_Yk_tau = tau;
  }
  if (m_Zk_tau != tau) {
    delete[] m_Zk;
    m_Zk = NULL;
    m_Zk_tau = tau;
  }
}

//-----------------------------------------------------------------------------

/// \brief  Initialises the pulse shape vector gk for the specified timing offset.
/// \param  tau       timing offset
///
void TimingEstimator::initgk(VALTYPE tau)
{
  int n;
  VALTYPE Tpshape;
  Tpshape = m_pulseshape->getduration();
  
  m_minn = ceil((-Tpshape/2-tau)/m_Delta)+1;
  m_maxn = floor((Tpshape/2-tau)/m_Delta)+1;
  
  if (m_gk_tau != tau) {
    delete[] m_gk;
    m_gk = NULL;
    // store tau for which m_gk will be computed
    m_gk_tau = tau;
  }
  
  if (!m_gk) {
    m_gk = new VALTYPE [(unsigned int)(m_maxn-m_minn+1)];
  }
  
  for (n=m_minn; n<=m_maxn; n++) {
    m_gk[n-m_minn] = m_pulseshape->pulseshape(-(n-1)*m_Delta - tau);
  }
}

//-----------------------------------------------------------------------------

/// \brief  Returns bk for the specified time index.
/// \param  k       time index
///
complexT TimingEstimator::getbk(int k)
{
  unsigned int l, minl, maxl, bkidx;
  int gkidxoffset;

  // allocate memory and initialise if called for the first time
  if (!m_bk) {
    m_bk = new complexT [m_K+m_a*m_q*(m_maxI-m_minI)];
    m_bkidx = new bool [m_K+m_a*m_q*(m_maxI-m_minI)];
    for (l=0; l < m_K+m_a*m_q*(m_maxI-m_minI); l++) {
      m_bk[l] = complexT (0,0);
      m_bkidx[l] = false; // boolean indicator if bk has been computed
    }
  }
  
  bkidx = k - m_a*m_q*m_minI - 1;
  
  if (!m_bkidx[bkidx]) {
    minl = MAX((int)  m_c-1,      k-m_maxn-1);
    maxl = MIN((int) (m_c*m_N)-1, k-m_minn-1);
    
    gkidxoffset = -m_minn-1; // map to indices starting with 0
        
    for (l=minl; l <= maxl; l++) {
      m_bk[bkidx] += m_zk[l] * m_gk[gkidxoffset + k-l];
    }
    
    m_bkidx[bkidx] = true;
  }
  
  return m_bk[bkidx];
}

//-----------------------------------------------------------------------------

/// \brief  Computes Yk for the specified timing offset and stores Yk in member
///         variable m_Yk.
/// \param  K         number of time grid points
/// \param  tau       timing offset
///
void TimingEstimator::getYk(unsigned int K, VALTYPE tau)
{
  int i;
  unsigned int k;
  int iMin, iMax;
    
  if (!m_Yk) {
    m_Yk = new complexT [K];
    for (k=0; k<K; k++) {
      m_Yk[k] = complexT (0,0);
    }
  }
  
  // Define these signed start/end indices here to work-around signed compare warning
  iMin = -m_a*m_q*m_maxP;
  iMax = -m_a*m_q*m_minP;

  // linear convolution
  for (k=1; k<=K; k++) {
    //for (i=-m_a*m_q*m_maxP; i<=-m_a*m_q*m_minP; i++) {
    for (i=iMin; i <= iMax; i++) {
      if (m_hP[i + m_a*m_q*m_maxP] != complexT (0,0)) {
        m_Yk[k-1] += m_hP[i + m_a*m_q*m_maxP] * getbk(k - i);
      }
    }
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
/// \brief  Returns Yk for the specified timing offset. Assumes a single time 
///         grid point.
/// \param  tau       timing offset
///
complexT TimingEstimator::getYk(VALTYPE tau)
{
  unsigned int K;
  K = 1;
  
  getYk(K, tau);
  return m_Yk[0];
}

//-----------------------------------------------------------------------------

/// \brief  Computes Zk for the specified timing offset and stores Zk in member
///         variable m_Zk.
/// \param  K         number of time grid points
/// \param  tau       timing offset
///
void TimingEstimator::getZk(unsigned int K, VALTYPE tau)
{
  unsigned int k, aqi;
  int i;
    
  if (!m_Zk) {
    m_Zk = new VALTYPE [K]();
  }
  
  if (m_maxD-m_minD+1 == m_numData) {    
    // contiguous data
    
    for (k=1; k<=MIN(m_a*m_q,K); k++) {
      for (aqi=m_a*m_q*m_minD; aqi<=m_maxD*m_a*m_q; aqi+=m_a*m_q) {
        m_Zk[k-1] += abs(getbk(k + aqi));
      }
    }
    // recursively compute Zk using sliding window
    for (k=m_a*m_q+1; k<=K; k++) {
      m_Zk[k-1] = m_Zk[k-1 - m_a*m_q] - abs(getbk(k + (m_minD-1)*m_a*m_q)) + abs(getbk(k + m_maxD*m_a*m_q));
    }
    
  } else {
    // non-contiguous data
    int iMin, iMax;
    // Define these signed start/end indices here to work-around signed compare warning
    iMin = -m_a*m_q*m_maxD;
    iMax = -m_a*m_q*m_minD;
    // linear convolution
    for (k=1; k<=K; k++) {
//      for (i=-m_a*m_q*m_maxD; i<=-m_a*m_q*m_minD; i++) {
      for (i=iMin; i <= iMax; i++) {
        if (m_hD[i + m_a*m_q*m_maxD]) {
          m_Zk[k-1] += abs(getbk(k - i));
          //m_Zk[k-1] += m_hD[i + m_a*m_q*m_maxD] * abs(getbk(k - i));
        }
      }
    }
    
  }  
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
/// \brief  Returns Zk for the specified timing offset. Assumes a single time 
///         grid point.
/// \param  tau       timing offset
///
VALTYPE TimingEstimator::getZk(VALTYPE tau)
{
  unsigned int K;
  K = 1;
    
  getZk(K, tau);
  return m_Zk[0];
}

//-----------------------------------------------------------------------------

/// \brief  Returns value of the sum of squares objective function for the 
///         specified timing offset.
/// \param  tau       timing offset
///
VALTYPE TimingEstimator::getSS(VALTYPE tau)
{
  checkVectorsUpToDate(tau);
  
  return -pow(getZk(tau) + abs(getYk(tau)),2)/m_L;
}

//-----------------------------------------------------------------------------

/// \brief  Refines the timing offset estimate given a coarse estimate. The 
///         algorithm uses Brent's golden section search and inverse parabolic 
///         interpolation algorithm. The default tolerance is 1e-4, but can be 
///         changed using setTolerance(tol).
/// \param  tau       coarse timing offset
///
void TimingEstimator::refine(VALTYPE tau)
{
  VALTYPE a, b;
  
  a = tau - m_Delta;
  b = tau + m_Delta;
  
  // run Brent's golden section search and inverse parabolic interpolation algorithm
  m_funcVal = local_min(a, b, m_refine_tol, m_timingEstimate);

  return;
}

//-----------------------------------------------------------------------------

/// \brief  Sets the tolerance for Brent's algorithm.
///
void TimingEstimator::setTolerance(VALTYPE tol)
{
  m_refine_tol = tol;
}

//-----------------------------------------------------------------------------

/// \brief  Returns the timing estimate.
///
VALTYPE TimingEstimator::getTimingEstimate()
{
  return m_timingEstimate;
}

//-----------------------------------------------------------------------------

/// \brief  Returns the value of the objective function at the estimated timing 
///         offset.
///
VALTYPE TimingEstimator::getFuncValue()
{
  return m_funcVal;
}

//-----------------------------------------------------------------------------

// BELOW FUNCTIONS COPIED FROM brent.cpp (GNU LGPL license). Modified to directly
// minimise getSS() rather than a generic function double (*f) (double) as C++ 
// member can't be passed as pointers. Also changed type of all variables from
// double to VALTYPE. The original code can be found at:
//
//      http://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.html

//****************************************************************************80

VALTYPE TimingEstimator::local_min ( VALTYPE a, VALTYPE b, VALTYPE t,
                  VALTYPE &x )

//****************************************************************************80
//
//  Purpose:
//
//    LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
//
//  Discussion:
//
//    The method used is a combination of golden section search and
//    successive parabolic interpolation.  Convergence is never much slower
//    than that for a Fibonacci search.  If F has a continuous second
//    derivative which is positive at the minimum (which is not at A or
//    B), then convergence is superlinear, and usually of the order of
//    about 1.324....
//
//    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
//    F is never evaluated at two points closer than TOL.
//
//    If F is a unimodal function and the computed values of F are always
//    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
//    LOCAL_MIN approximates the abscissa of the global minimum of F on the
//    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
//
//    If F is not unimodal, then LOCAL_MIN may approximate a local, but
//    perhaps non-global, minimum to the same accuracy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2011
//
//  Author:
//
//    Original FORTRAN77 version by Richard Brent.
//    C++ version by John Burkardt.
//    Modifications by John Denker.
//
//  Reference:
//
//    Richard Brent,
//    Algorithms for Minimization Without Derivatives,
//    Dover, 2002,
//    ISBN: 0-486-41998-3,
//    LC: QA402.5.B74.
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, double T, a positive absolute error tolerance.
//
//    Input, func_base& F, a user-supplied c++ functor whose
//    local minimum is being sought.  The input and output
//    of F() are of type double.
//
//    Output, double &X, the estimated value of an abscissa
//    for which F attains a local minimum value in [A,B].
//
//    Output, double LOCAL_MIN, the value F(X).
//
{
  VALTYPE c;
  VALTYPE d;
  VALTYPE e;
  VALTYPE eps;
  VALTYPE fu;
  VALTYPE fv;
  VALTYPE fw;
  VALTYPE fx;
  VALTYPE m;
  VALTYPE p;
  VALTYPE q;
  VALTYPE r;
  VALTYPE sa;
  VALTYPE sb;
  VALTYPE t2;
  VALTYPE tol;
  VALTYPE u;
  VALTYPE v;
  VALTYPE w;
  //
  //  C is the square of the inverse of the golden ratio.
  //
  c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );
  
  eps = sqrt ( r8_epsilon ( ) );
  
  sa = a;
  sb = b;
  x = sa + c * ( b - a );
  w = x;
  v = w;
  e = 0.0;
  fx = getSS(x);
  fw = fx;
  fv = fw;
  
  for ( ; ; )
  {
    m = 0.5 * ( sa + sb ) ;
    tol = eps * r8_abs ( x ) + t;
    t2 = 2.0 * tol;
    //
    //  Check the stopping criterion.
    //
    if ( r8_abs ( x - m ) <= t2 - 0.5 * ( sb - sa ) )
    {
      break;
    }
    //
    //  Fit a parabola.
    //
    r = 0.0;
    q = r;
    p = q;
    
    if ( tol < r8_abs ( e ) )
    {
      r = ( x - w ) * ( fx - fv );
      q = ( x - v ) * ( fx - fw );
      p = ( x - v ) * q - ( x - w ) * r;
      q = 2.0 * ( q - r );
      if ( 0.0 < q )
      {
        p = - p;
      }
      q = r8_abs ( q );
      r = e;
      e = d;
    }
    
    if ( r8_abs ( p ) < r8_abs ( 0.5 * q * r ) &&
        q * ( sa - x ) < p &&
        p < q * ( sb - x ) )
    {
      //
      //  Take the parabolic interpolation step.
      //
      d = p / q;
      u = x + d;
      //
      //  F must not be evaluated too close to A or B.
      //
      if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
      {
        if ( x < m )
        {
          d = tol;
        }
        else
        {
          d = - tol;
        }
      }
    }
    //
    //  A golden-section step.
    //
    else
    {
      if ( x < m )
      {
        e = sb - x;
      }
      else
      {
        e = sa - x;
      }
      d = c * e;
    }
    //
    //  F must not be evaluated too close to X.
    //
    if ( tol <= r8_abs ( d ) )
    {
      u = x + d;
    }
    else if ( 0.0 < d )
    {
      u = x + tol;
    }
    else
    {
      u = x - tol;
    }
    
    fu = getSS( u );
    //
    //  Update A, B, V, W, and X.
    //
    if ( fu <= fx )
    {
      if ( u < x )
      {
        sb = x;
      }
      else
      {
        sa = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
    else
    {
      if ( u < x )
      {
        sa = u;
      }
      else
      {
        sb = u;
      }
      
      if ( fu <= fw || w == x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu <= fv || v == x || v== w )
      {
        v = u;
        fv = fu;
      }
    }
  }
  return fx;
}
//****************************************************************************80

VALTYPE TimingEstimator::r8_abs ( VALTYPE x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, VALTYPE X, the quantity whose absolute value is desired.
//
//    Output, VALTYPE R8_ABS, the absolute value of X.
//
{
  VALTYPE value;
  
  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

VALTYPE TimingEstimator::r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 round off unit.
//
//  Discussion:
//
//    R8_EPSILON is a number R which is a power of 2 with the property that,
//    to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, VALTYPE R8_EPSILON, the VALTYPE precision round-off unit.
//
{
  VALTYPE r;
  
  r = 1.0;
  
  while ( 1.0 < ( VALTYPE ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }
  
  return ( 2.0 * r );
}
//****************************************************************************80

