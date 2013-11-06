//
//  QAMModulatorV2.cpp
//  QAMModulator
//
//  Created by Gottfried Lechner on 16/03/12.
//
//  Converted to use vector and complex templates and
//  use ARM NEON intrinsics for the case of 4 symbols and floats
//  Gottfried Lechner 09/02/13
//
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "QAMModulatorV2.h"

#ifdef __ARM_NEON__
  #include "arm_neon.h"
  #include "neon_mathfun.h"
#endif

#define   POW2(x) ((x)*(x))

#ifdef FAST_MATH
  #include  "fastexp.h"
  #include  "fastlog.h"
        
  #define   LOG(x)  fasterlog(x)
  #define   EXP(x)  fasterexp(x)
#else
  #define   LOG(x)  log(x)
  #define   EXP(x)  exp(x)
#endif

#ifdef _WIN32
  #define log2(x) 	(log(x)/log((MSGTYPE)2.0))
#endif


/// \brief	defines the maximum absolute value of L-values used
///
#define MAXLLR	100

/// \brief	returns minimum of a and b
///
#define MIN(a,b)			(((a) < (b)) ? (a) : (b))

/// \brief	returns maximum of a and b
///
#define MAX(a,b)			(((a) > (b)) ? (a) : (b))

/// \brief	clips a at magnitude b
///
#define CLIP(a,b)			MAX(MIN(a,b), -b)

// macros for memory access
#define PCH(time,symbol)            (m_pch[(time)*m_numSymbols + (symbol)])
#define LA(time,bit)                (m_La[(time)*m_bitsperSymbol + (bit)])

/// \brief  Standard constructor.
///
CQAMModulator::CQAMModulator(void)
{ 
  m_bitsperSymbol = 0;
  m_numSections   = 0;
  m_numSymbols    = 0;
  m_apriset       = 0;
}

/// \brief  Standard destructor.
///
CQAMModulator::~CQAMModulator(void)
{
}

/// \brief  Initialises modulator.
/// \param  numSections   length of received vector
/// \param  numSymbols    number of symbols per time slot
/// \param  SymReal       vector of real parts of constellation
/// \param  SymImag       vector of imaginary parts of constellation
///
/// Copies the constellation definition to the local memory, and calls
/// the memory allocation.
void CQAMModulator::init(unsigned int numSections,
                         unsigned int numSymbols,
                         std::vector<complex> Sym)
{   
  m_numSections   = numSections;
  m_numSymbols    = numSymbols;
  m_bitsperSymbol = log2((MSGTYPE)numSymbols);
  m_apriset       = 0;
  
  m_Sym           = Sym;
  
  m_pch.resize(m_numSections*m_numSymbols);
}

/// \brief  Sets the symbol probabilities based on the received values.
/// \param  yReal   vector of real parts of received signal
/// \param  yImag   vector of real parts of received signal
/// \param  sigma2  vector of noise variances
void CQAMModulator::set_ch(std::vector<complex> y, std::vector<MSGTYPE> sigma2)
{
  unsigned int  time, symbol;
  
  for(time=0;time<m_numSections;time++)
    for(symbol=0;symbol<m_numSymbols;symbol++)
      PCH(time,symbol) = EXP(-(POW2(m_Sym[symbol].real()-y[time].real()) + POW2(m_Sym[symbol].imag()-y[time].imag()))/sigma2[time]);
}

/// \brief  Sets the symbol probabilities based on the received values with uniform noise.
/// \param  yReal   vector of real parts of received signal
/// \param  yImag   vector of real parts of received signal
/// \param  sigma2  noise variance
void CQAMModulator::set_ch(std::vector<complex> y, MSGTYPE sigma2)
{
  unsigned int  time, symbol;
  MSGTYPE       sigma2inv;
  
  sigma2inv = (MSGTYPE)(1.0)/sigma2;
  
  
#ifdef __ARM_NEON__
  if((m_numSymbols!=4) || (sizeof(MSGTYPE)!=4))
  {
  	std::cerr << "NEON only supported for 4 symbols and floats but " << m_numSymbols << " symbols and " << sizeof(MSGTYPE) << " byte messages are configured.";
  	exit(0);
  }
  
  // load all symbols with real and imaginary part in two variables
  float32x4_t	symreal, symimag, yreal, yimag, dreal, dimag, dist, prob;
  for(symbol=0;symbol<4;symbol++)
  {
	symreal = vld1q_lane_f32(&(m_Sym[symbol].real()), symreal, symbol);
	symimag = vld1q_lane_f32(&(m_Sym[symbol].imag()), symimag, symbol);
  }

  for(time=0;time<m_numSections;time++)
  {
    // load symbol
    yreal = vld1q_dup_f32(&(y[time].real()));
    yimag = vld1q_dup_f32(&(y[time].imag()));
    
    dreal = vsubq_f32(symreal, yreal);
    dimag = vsubq_f32(symimag, yimag);
    
    dist  = vmulq_f32(dreal, dreal);
    dist  = vmlaq_f32(dist, dimag, dimag);
    
    dist  = vmulq_n_f32(dist, -sigma2inv);
    
    prob  = exp_ps(dist);
    
    vst1q_f32(&(PCH(time,0)), prob);
  }	
#else  
  for(time=0;time<m_numSections;time++)
    for(symbol=0;symbol<m_numSymbols;symbol++)
      PCH(time,symbol) = EXP(-(POW2(m_Sym[symbol].real()-y[time].real()) + POW2(m_Sym[symbol].imag()-y[time].imag()))*sigma2inv);
#endif    
}

/// \brief  Sets the a-priori LLRs.
/// \param  La   a-priori LLRs (dim: bitspersymbol x numsections)
void CQAMModulator::set_apri(std::vector<MSGTYPE> La)
{
  m_La = La;
}

/// \brief  Modulates binary data.
/// \param  data   binary input vector
/// \param  yReal  vector of real parts of modulated signal
/// \param  yImag  vector of real parts of modulated signal
void CQAMModulator::modulate(unsigned int *data, std::vector<complex>& y)
{
  unsigned int  time, symbol, bit;
  
  for(time=0;time<m_numSections;time++)
  {  
    symbol = 0;
    for(bit=0;bit<m_bitsperSymbol;bit++)
      symbol = (symbol<<1) | data[time*m_bitsperSymbol + bit];
  
    y[time] = m_Sym[symbol];
  }
}

/// \brief  Demodulate and returns LLRs for input symbols.
/// \param  Lout   output LLRs (dim: bitspersymbol x numsections)
void CQAMModulator::demodulate(std::vector<MSGTYPE>& Lout)
{
  unsigned int  time, symbol, bit, bita, b;
  MSGTYPE       Ps, p[2];
  
  #ifdef __ARM_NEON__
    float32x4_t	x, xmin, xmax;
    float32_t 	maxllr = MAXLLR;
    xmax   = vld1q_dup_f32(&maxllr);
    maxllr = -maxllr;
    xmin   = vld1q_dup_f32(&maxllr);
  #endif  
  
  // loop through received vector
  for(time=0;time<m_numSections;time++)
  {
    // compute extrinsic for every bit
    for(bit=0;bit<m_bitsperSymbol;bit++)
    {
      p[0] = 0.0;
      p[1] = 0.0;
      
      // sum over symbols
      for(symbol=0;symbol<m_numSymbols;symbol++)
      {
        if(m_apriset)
        {        
          // start with channel probability
          // we already computed the EXP in set_chan which has to be undone in case
          // a-priori information is used
          Ps = LOG(PCH(time,symbol));
          // include a-priori information of all other bits
          for(bita=0;bita<m_bitsperSymbol;bita++)
          {
            if(bita != bit)
            {
              b   = (symbol & (1<<bita))>>bita;
              Ps += (1-2*(MSGTYPE)(b))*LA(time,m_bitsperSymbol-bita-1)/2 - LOG(EXP(LA(time,m_bitsperSymbol-bita-1)/2)+EXP(-LA(time,m_bitsperSymbol-bita-1)/2));
            }
          }
          p[(symbol & (1<<bit))>>bit] += EXP(Ps);
        }
        else
          p[(symbol & (1<<bit))>>bit] += PCH(time,symbol);
      }

#ifdef __ARM_NEON__
    Lout[(time)*m_bitsperSymbol + (m_bitsperSymbol-bit-1)] = p[0]/p[1];
#else
    Lout[(time)*m_bitsperSymbol + (m_bitsperSymbol-bit-1)] = CLIP(LOG(p[0]/p[1]),MAXLLR);
#endif      
    }
  }
  
#ifdef __ARM_NEON__  
  for(time=0;time<m_numSections;time+=2)
  {
    x = vld1q_f32(&(Lout[time*m_bitsperSymbol]));
    x = log_ps(x);
    x = vmaxq_f32(x, xmin);
    x = vminq_f32(x, xmax);
    vst1q_f32(&(Lout[time*m_bitsperSymbol]), x);
  }  
#endif
}


/// \brief  Remodulates and returns the soft output symbols.
/// \param  yReal   vector of real parts of remodulated signal
/// \param  yImag   vector of real parts of remodulated signal
void CQAMModulator::remodulate(std::vector<complex>& y, std::vector<MSGTYPE>& ResVar)
{
  unsigned int  time, symbol, bit, b;
  MSGTYPE       *Ps, sum;
  
  Ps = (MSGTYPE*)malloc(m_numSymbols * sizeof(MSGTYPE));
  
  // loop through received vector
  for(time=0;time<m_numSections;time++)
  {   
    y[time] = (0.0, 0.0);
    
    // compute probability for every symbol
    sum = 0;
    for(symbol=0;symbol<m_numSymbols;symbol++)
    {
      Ps[symbol] = 0;
        
      // include a-priori information of all bits
      for(bit=0;bit<m_bitsperSymbol;bit++)
      {
        b           = (symbol & (1<<bit))>>bit;
        Ps[symbol] += (1-2*(MSGTYPE)(b))*LA(time,m_bitsperSymbol-bit-1)/2 - LOG(EXP(LA(time,m_bitsperSymbol-bit-1)/2)+EXP(-LA(time,m_bitsperSymbol-bit-1)/2));
      }
      
      // convert log-prob to probabilities and sum them up (used for normalisation later)
      Ps[symbol] = EXP(Ps[symbol]);
      sum       += Ps[symbol];
    }

    // normalise
    for(symbol=0;symbol<m_numSymbols;symbol++)
      Ps[symbol] /= sum;    
    
    // add contribution to remodulated signal
    for(symbol=0;symbol<m_numSymbols;symbol++)
      y[time] += Ps[symbol]*m_Sym[symbol];

    // compute the variance of the residual signal
    ResVar[time] = 0;
    for(symbol=0;symbol<m_numSymbols;symbol++)
      ResVar[time] += Ps[symbol] * (POW2(y[time].real()-m_Sym[symbol].real()) + POW2(y[time].imag()-m_Sym[symbol].imag()));
  }
  
  free(Ps);
}

/// \brief  Returns the length of the received signal.
/// 
unsigned int CQAMModulator::getlength()
{
  return m_numSections;
}

/// \brief  Returns the number of bits per symbol.
/// 
unsigned int CQAMModulator::getbitspersymbol()
{
  return m_bitsperSymbol;
}
