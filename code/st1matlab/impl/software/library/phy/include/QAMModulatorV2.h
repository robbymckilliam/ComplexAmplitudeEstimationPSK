//
//  QAMModulatorV2.h
//  QAMModulator
//
//  Created by Gottfried Lechner on 09/02/13.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#ifndef QAMModulator_QAMModulator_h
#define QAMModulator_QAMModulator_h

#ifndef MSGTYPE
  #define MSGTYPE double
#endif
      
#include <complex>
#include <vector>
typedef std::complex<MSGTYPE> complex;
        
class CQAMModulator
{

public:
  CQAMModulator(void);
  ~CQAMModulator(void);

  void          init(unsigned int numSections,
                     unsigned int numSymbols,
                     std::vector<complex> Sym);
    
  void          set_ch(std::vector<complex> y, std::vector<MSGTYPE> sigma2);
  void          set_ch(std::vector<complex> y, MSGTYPE sigma2);
  
  void          set_apri(std::vector<MSGTYPE> La);
  
  void          modulate(unsigned int *data, std::vector<complex>& y);
  void          demodulate(std::vector<MSGTYPE>& Lout);
  void          remodulate(std::vector<complex>& y, std::vector<MSGTYPE>& ResVar);
  
  unsigned int  getlength();
  unsigned int  getbitspersymbol();
  
private:
  unsigned int  m_numSections;
  unsigned int  m_numSymbols;
  unsigned int  m_bitsperSymbol;
  unsigned int  m_apriset;
  
  std::vector<complex>  m_Sym;
  std::vector<MSGTYPE>  m_La;
  std::vector<MSGTYPE>  m_pch;
};

#endif
