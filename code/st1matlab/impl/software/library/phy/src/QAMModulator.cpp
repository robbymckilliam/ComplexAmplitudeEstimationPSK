//
//  QAMModulator.cpp
//  QAMModulator
//
//  Created by Gottfried Lechner on 16/03/12.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#include <iostream>//iostream.h
#include <math.h>
#include <stdlib.h>
#include "QAMModulator.h"

#ifdef FAST_MATH
#include "fastexp.h"
#include "fastlog.h"
#endif

#ifdef _WIN32
#define log2(x) 			(log(x)/log((MSGTYPE)2.0))
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
#define PA(time,symbol)             (m_pa[(time)*m_numSymbols + (symbol)])

#define LA(time,bit)                (m_La[(time)*m_bitsperSymbol + (bit)])
#define LOUT(time,bit)              (Lout[(time)*m_bitsperSymbol + (bit)])

/// \brief  Standard constructor.
///
CQAMModulator::CQAMModulator(void)
{
  m_SymReal       = NULL;
  m_SymImag       = NULL;
  m_pch           = NULL;
  m_La            = NULL;
  
  m_bitsperSymbol = 0;
  m_numSections   = 0;
  m_numSymbols    = 0;
}

/// \brief  Standard destructor.
///
CQAMModulator::~CQAMModulator(void)
{
  freeMem();

  m_bitsperSymbol = 0;
  m_numSections   = 0;
  m_numSymbols    = 0;
}

/// \brief  Memory allocation.
///
/// Allocates memory to hold the branch metrics, state probabilities and the
/// trellis structure.
void CQAMModulator::allocateMem()
{
  unsigned int  n;
  
  m_SymReal = (MSGTYPE*)malloc(m_numSymbols * sizeof(MSGTYPE));
  m_SymImag = (MSGTYPE*)malloc(m_numSymbols * sizeof(MSGTYPE));
  m_pch     = (MSGTYPE*)malloc(m_numSections * m_numSymbols * sizeof(MSGTYPE));
  m_La      = (MSGTYPE*)malloc(m_numSections * m_bitsperSymbol * sizeof(MSGTYPE));
  
  // Initialise with zero a-priori values (in case the class is just used as
  // a demodulator).
  for(n=0;n<m_numSections*m_bitsperSymbol;n++)
    m_La[n] = 0.0;
}

/// \brief  Frees allocated memory.
///
void CQAMModulator::freeMem()
{
  free(m_SymReal);
  free(m_SymImag);
  free(m_pch);
  free(m_La);
  
  m_SymReal = NULL;
  m_SymImag = NULL;
  m_pch     = NULL;
  m_La      = NULL;  
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
                         MSGTYPE      *SymReal,
                         MSGTYPE      *SymImag)
{
  unsigned int  symbol;
  
  freeMem();
  
  m_numSections   = numSections;
  m_numSymbols    = numSymbols;
  m_bitsperSymbol = log2((MSGTYPE)numSymbols);
  
  allocateMem();
  
  for(symbol=0;symbol<m_numSymbols;symbol++)
  {
    m_SymReal[symbol] = SymReal[symbol];
    m_SymImag[symbol] = SymImag[symbol];
  }
}

/// \brief  Sets the symbol probabilities based on the received values.
/// \param  yReal   vector of real parts of received signal
/// \param  yImag   vector of real parts of received signal
void CQAMModulator::set_ch(MSGTYPE *yReal, MSGTYPE *yImag, MSGTYPE *sigma2)
{
  unsigned int  time, symbol;
  
  for(time=0;time<m_numSections;time++)
    for(symbol=0;symbol<m_numSymbols;symbol++)
#ifdef FAST_MATH
	PCH(time,symbol) = -((m_SymReal[symbol]-yReal[time])*(m_SymReal[symbol]-yReal[time]) + (m_SymImag[symbol]-yImag[time])*(m_SymImag[symbol]-yImag[time]))/sigma2[time];
#else
	PCH(time,symbol) = -(pow(m_SymReal[symbol]-yReal[time], 2.0) + pow(m_SymImag[symbol]-yImag[time], 2.0))/sigma2[time];
#endif		
}

/// \brief  Sets the a-priori LLRs.
/// \param  La   a-priori LLRs (dim: bitspersymbol x numsections)
void CQAMModulator::set_apri(MSGTYPE *La)
{
  unsigned int  n;
  
  for(n=0;n<m_bitsperSymbol*m_numSections;n++)
    m_La[n] = La[n];
}

/// \brief  Modulates binary data.
/// \param  data   binary input vector
/// \param  yReal  vector of real parts of modulated signal
/// \param  yImag  vector of real parts of modulated signal
void CQAMModulator::modulate(unsigned int *data, MSGTYPE *yReal, MSGTYPE *yImag)
{
  unsigned int  time, symbol, bit;
  
  for(time=0;time<m_numSections;time++)
  {  
    symbol = 0;
    for(bit=0;bit<m_bitsperSymbol;bit++)
      symbol = (symbol<<1) | data[time*m_bitsperSymbol + bit];
    
    yReal[time] = m_SymReal[symbol];
    yImag[time] = m_SymImag[symbol];    
  }
}

/// \brief  Demodulate and returns LLRs for input symbols.
/// \param  Lout   output LLRs (dim: bitspersymbol x numsections)
void CQAMModulator::demodulate(MSGTYPE *Lout)
{
  unsigned int  time, symbol, bit, bita, b;
  MSGTYPE       Ps, p[2];
  
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
        // start with channel probability
        Ps = PCH(time,symbol);
        
        // include a-priori information of all other bits
        for(bita=0;bita<m_bitsperSymbol;bita++)
        {
          if(bita != bit)
          {
            b   = (symbol & (1<<bita))>>bita;
#ifdef FAST_MATH
	    Ps += (1-2*(MSGTYPE)(b))*LA(time,m_bitsperSymbol-bita-1)/2 - fasterlog(fasterexp(LA(time,m_bitsperSymbol-bita-1)/2)+fasterexp(-LA(time,m_bitsperSymbol-bita-1)/2));
#else
	    Ps += (1-2*(MSGTYPE)(b))*LA(time,m_bitsperSymbol-bita-1)/2 - log(exp(LA(time,m_bitsperSymbol-bita-1)/2)+exp(-LA(time,m_bitsperSymbol-bita-1)/2));
#endif
          }
        }
#ifdef FAST_MATH
	p[(symbol & (1<<bit))>>bit] += fasterexp(Ps);
#else
        p[(symbol & (1<<bit))>>bit] += exp(Ps);
#endif
      }
#ifdef FAST_MATH
      LOUT(time,m_bitsperSymbol-bit-1) = CLIP(fasterlog(p[0]/p[1]),MAXLLR);
#else
      LOUT(time,m_bitsperSymbol-bit-1) = CLIP(log(p[0]/p[1]),MAXLLR);
#endif
    }
  }
}

/// \brief  Remodulates and returns the soft output symbols.
/// \param  yReal   vector of real parts of remodulated signal
/// \param  yImag   vector of real parts of remodulated signal
void CQAMModulator::remodulate(MSGTYPE *yReal, MSGTYPE *yImag, MSGTYPE *ResVar=NULL)
{
  unsigned int  time, symbol, bit, b;
  MSGTYPE       *Ps, sum;
  
  Ps = (MSGTYPE*)malloc(m_numSymbols * sizeof(MSGTYPE));
  
  // loop through received vector
  for(time=0;time<m_numSections;time++)
  {   
    yReal[time] = 0.0;
    yImag[time] = 0.0;
    
    // compute probability for every symbol
    sum = 0;
    for(symbol=0;symbol<m_numSymbols;symbol++)
    {
      Ps[symbol] = 0;
        
      // include a-priori information of all bits
      for(bit=0;bit<m_bitsperSymbol;bit++)
      {
        b           = (symbol & (1<<bit))>>bit;
        Ps[symbol] += (1-2*(MSGTYPE)(b))*LA(time,m_bitsperSymbol-bit-1)/2 - log(exp(LA(time,m_bitsperSymbol-bit-1)/2)+exp(-LA(time,m_bitsperSymbol-bit-1)/2));
      }
      
      // convert log-prob to probabilities and sum them up (used for normalisation later)
      Ps[symbol] = exp(Ps[symbol]);
      sum += Ps[symbol];
    }

    // normalise
    for(symbol=0;symbol<m_numSymbols;symbol++)
      Ps[symbol] /= sum;    
    
    for(symbol=0;symbol<m_numSymbols;symbol++)
    {
      // add contribution to remodulated signal
      yReal[time] += Ps[symbol]*m_SymReal[symbol];
      yImag[time] += Ps[symbol]*m_SymImag[symbol];
    }
    
    if(ResVar)
    {
      // compute the variance of the residual signal
      ResVar[time] = 0;
      for(symbol=0;symbol<m_numSymbols;symbol++)
        ResVar[time] += Ps[symbol] * (pow(yReal[time]-m_SymReal[symbol],2.0) + pow(yImag[time]-m_SymImag[symbol],2.0));
    }
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
