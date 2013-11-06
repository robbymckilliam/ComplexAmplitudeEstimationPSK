//
//  TimingEstimator.h
//  TimingEstimator
//
//  Created by Andre Pollok on 31/05/2012.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#ifndef TimingEstimator_TimingEstimator_h
#define TimingEstimator_TimingEstimator_h

#ifndef VALTYPE
#define VALTYPE double
#endif
#ifndef MSGTYPE
#define MSGTYPE VALTYPE
#endif

typedef std::complex<VALTYPE> complexT;

class PulseShape
{
public:
  PulseShape(VALTYPE Tsymb, VALTYPE OS, VALTYPE duration);
  ~PulseShape(void) {};
  
  VALTYPE getTsymb();
  VALTYPE getOS();
  VALTYPE getduration();
  virtual VALTYPE pulseshape(VALTYPE time) = 0;
  
protected:
  VALTYPE m_Tsymb;
  VALTYPE m_OS;
  VALTYPE m_duration;
};

class RootRaisedCosine : public PulseShape
{
public:
  RootRaisedCosine(VALTYPE rolloff, VALTYPE Tsymb, VALTYPE OS, VALTYPE duration);  
  ~RootRaisedCosine(void);
  
//  VALTYPE getduration();
  
  VALTYPE pulseshape(VALTYPE time);
  
protected:
  VALTYPE m_rolloff;
  
  VALTYPE m_fc, m_fs, m_tol;
};

class TimingEstimator
{

public:
  TimingEstimator(VALTYPE          taumin,
                  VALTYPE          taumax,
                  PulseShape*      pulseshape,
                  unsigned int*    D,
                  unsigned int     numData,
                  complexT*        pilot,
                  unsigned int*    P,
                  unsigned int     numPilots,
                  unsigned int     a=4);
  
  ~TimingEstimator(void);
  
  void             runEstimator(complexT* rxSamp, unsigned int numSamp);
  
  VALTYPE          getTimingEstimate();
  VALTYPE          getFuncValue();
  void             setTolerance(VALTYPE);
  
// ONLY FOR DEBUGGING //
  inline VALTYPE*  getgk() { return m_gk; };
  inline complexT* getzk() { return m_zk; };
  inline complexT* getbkpub() { return m_bk; };
  inline complexT* getYk() { return m_Yk; };
  inline VALTYPE*  getZk() { return m_Zk; };
// END --- ONLY FOR DEBUGGING //
  
private:  
  unsigned int     m_p;
  unsigned int     m_q;
  VALTYPE          m_T;
  VALTYPE          m_Ts;
  PulseShape*      m_pulseshape;
  unsigned int     m_c;
  unsigned int     m_a;
  VALTYPE          m_Delta;
  VALTYPE*         m_taugrid;
  VALTYPE          m_gk_tau; // last tau for which gk, bk, Yk and Zk have been computed
  VALTYPE          m_bk_tau;
  VALTYPE          m_Yk_tau;
  VALTYPE          m_Zk_tau;

  unsigned int     m_K;
  unsigned int     m_L;
  unsigned int     m_minD;
  unsigned int     m_maxD;
  unsigned int     m_minP;
  unsigned int     m_maxP;
  unsigned int     m_minI;
  unsigned int     m_maxI;
  unsigned int     m_N;
  int              m_minn;
  int              m_maxn;

  unsigned int     m_numData;
  unsigned int*    m_D;
  bool*            m_hD;
  unsigned int     m_numPilots;
  unsigned int*    m_P;
  complexT*        m_pilot;
  complexT*        m_hP;
  
  complexT*        m_zk;
  complexT*        m_bk;
  bool*            m_bkidx;
  VALTYPE*         m_gk;
  VALTYPE*         m_Zk;
  complexT*        m_Yk;
  
  VALTYPE          m_refine_tol;
  
  VALTYPE          m_timingEstimate;
  VALTYPE          m_funcVal;
  
  void             checkVectorsUpToDate(VALTYPE tau);
  void             initzk(complexT*);
  void             initgk(VALTYPE);
  complexT         getbk(int);
  void             getYk(unsigned int, VALTYPE);
  complexT         getYk(VALTYPE);
  void             getZk(unsigned int, VALTYPE);
  VALTYPE          getZk(VALTYPE);
  VALTYPE          getSS(VALTYPE);
  void             refine(VALTYPE);
  
  // Brent's functions
  VALTYPE          local_min(VALTYPE, VALTYPE, VALTYPE,
                             VALTYPE&);
  VALTYPE          r8_abs(VALTYPE);
  VALTYPE          r8_epsilon();
};

#endif
