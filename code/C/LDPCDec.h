#ifndef LDPCDEC_H
#define LDPCDEC_H

#define DECODER_SPA		0
#define DECODER_MSA		1

#ifndef MSGTYPE
  #define MSGTYPE double
#endif

/// \brief		Main class for LDPC decoding
/// \author		Gottfried Lechner (gottfried.lechner@unisa.edu.au)
/// \version	3.0
/// \date     Created      : March 2008
/// \date     Last modified: April 2012
///
/// This class holds the LDPC decoder.
class CLDPCDec
{
public:
	CLDPCDec(void);
	~CLDPCDec(void);

	int						readDecoder(const char* filename);
	
  unsigned int  encodeLDGM(unsigned int *info, unsigned int *parity);
  unsigned int  encodeRA(unsigned int *info, unsigned int *parity);
	unsigned int	decode(MSGTYPE *Lch, MSGTYPE *Lapp, unsigned int maxit);

	void					set_param_clearmsg(bool clearmsg);
	void					set_param_method(unsigned int method);
  void          set_param_corrvec(unsigned int len, MSGTYPE *corrvec);

	MSGTYPE				syndromeInformation();	
  
	bool					dec_loaded();

	unsigned int	getN();
	unsigned int	getM();
	unsigned int	getK();
	unsigned int	getE();

private:
	unsigned int	m_N;
	unsigned int	m_M;
	unsigned int	m_E;
	unsigned int	m_g;
	unsigned int	*m_vardegree;
	unsigned int	*m_chkdegree;
	unsigned int	*m_interleaver;
	MSGTYPE				*m_c2v;
	MSGTYPE				*m_v2c;
	int						*m_hard;
  MSGTYPE				*m_Lsyn;

	bool					m_clearmsg;
	unsigned int	m_method;
	
	bool					m_decloaded;

  unsigned int  m_corrveclen;
  MSGTYPE       *m_corrvec;
  
	unsigned int	decodeSPA(MSGTYPE *Lch, MSGTYPE *Lapp, unsigned int maxit);
	unsigned int	decodeMSA(MSGTYPE *Lch, MSGTYPE *Lapp, unsigned int maxit);

	MSGTYPE				L2mutual(MSGTYPE *L, long N);  
  
	void					allocateMem();
	void					freeMem();
};

#endif
