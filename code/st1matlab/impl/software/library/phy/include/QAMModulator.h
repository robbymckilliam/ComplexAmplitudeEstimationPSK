//
//  QAMModulator.h
//  QAMModulator
//
//  Created by Gottfried Lechner on 16/03/12.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#ifndef QAMModulator_QAMModulator_h
#define QAMModulator_QAMModulator_h

#ifndef MSGTYPE
#define MSGTYPE double
#endif

class CQAMModulator
{

public:
  CQAMModulator(void);
  ~CQAMModulator(void);

  void          init(unsigned int numSections,
                     unsigned int numSymbols,
                     MSGTYPE      *SymReal,
                     MSGTYPE      *SymImag);
    
  void          set_ch(MSGTYPE *yReal, MSGTYPE *yImag, MSGTYPE *sigma2);
  void          set_apri(MSGTYPE *La);
  
  void          modulate(unsigned int *data, MSGTYPE *yReal, MSGTYPE *yImag);
  void          demodulate(MSGTYPE *Lout);
  void          remodulate(MSGTYPE *yReal, MSGTYPE *yImag, MSGTYPE *ResVar);
  
  unsigned int  getlength();
  unsigned int  getbitspersymbol();
    
private:
  unsigned int  m_numSections;
  unsigned int  m_numSymbols;
  unsigned int  m_bitsperSymbol;
  MSGTYPE        *m_SymReal, *m_SymImag;
  
  MSGTYPE       *m_pch;
  MSGTYPE       *m_pa;
  MSGTYPE       *m_La;
  
  void          allocateMem();
  void          freeMem();
};

#endif
