#include "LDPCDecV2.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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
  #define log2(x) 			(log(x)/log((MSGTYPE)2))
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

/// \brief	return sign of a
///
#define SIGN(a)				(((a)<0) ? (-1) : (1))

/// \brief	return hard decision of a
///
#define HD(a)				  (signbit(a))

//#define FABS(a)				(((a)<0) ? (-a) : (a))
#define FABS(a)				copysignf(a,1.0f)

/// \brief	Standard constructor.
///
CLDPCDec::CLDPCDec(void)
{
	m_N						= 0;
	m_M						= 0;
	m_E						= 0;
	m_g						= 0;
	m_vardegree		= NULL;
	m_chkdegree		= NULL;
	m_interleaver = NULL;
	m_v2c					= NULL;
	m_c2v					= NULL;
  m_tmp         = NULL;
	m_hard				= NULL;
  
  m_mf          = NULL;
  m_mb          = NULL;
  m_sf          = NULL;
  m_sb          = NULL;
  m_dcmax       = 0;

	m_clearmsg		= true;
	m_method			= DECODER_SPA;
	
	m_decloaded		= false;
  
  m_corrveclen  = 0;
  m_corrvec     = NULL;
}

/// \brief	Standard destructor.
///
CLDPCDec::~CLDPCDec(void)
{
	freeMem();
}

/// \brief  Allocate memory.
///
void CLDPCDec::allocateMem()
{
	freeMem();
	m_vardegree		= (unsigned int*)malloc(m_N * sizeof(unsigned int));
	m_chkdegree		= (unsigned int*)malloc(m_M * sizeof(unsigned int));
	m_interleaver	= (unsigned int*)malloc(m_E * sizeof(unsigned int));
	m_v2c					= (MSGTYPE*)malloc(m_E * sizeof(MSGTYPE));
	m_c2v					= (MSGTYPE*)malloc(m_E * sizeof(MSGTYPE));
 	m_tmp					= (MSGTYPE*)malloc(m_E * sizeof(MSGTYPE));
  // m_hard is used in the method encodeRA where it needs to be one element
  // larger than the number of edges. This is because the last edge of the
  // accumulator is ommited in the code definition.
	m_hard				= (int*)malloc((m_E+1) * sizeof(int));
	m_decloaded		= true;
}

/// \brief  Free memory.
///
void CLDPCDec::freeMem()
{
	free(m_vardegree);
	free(m_chkdegree);
	free(m_interleaver);
	free(m_v2c);
	free(m_c2v);
  free(m_tmp);
	free(m_hard);
  free(m_corrvec);
  free(m_mf);
  free(m_mb);
  free(m_sf);
  free(m_sb);
	
	m_vardegree		= NULL;
	m_chkdegree		= NULL;
	m_interleaver = NULL;
	m_v2c					= NULL;
	m_c2v					= NULL;
  m_tmp         = NULL;
	m_hard				= NULL;
	m_decloaded		= false;
  m_corrvec     = NULL;
  m_corrveclen  = 0;
  m_mf          = NULL;
  m_mb          = NULL;
  m_sf          = NULL;
  m_sb          = NULL;
}

/// \brief	Returns the length of the codeword.
///
unsigned int CLDPCDec::getN()
{
	return m_N;	
}

/// \brief	Returns the number of parity-checks.
///
unsigned int CLDPCDec::getM()
{
	return m_M;	
}

/// \brief	Returns the number of information bits.
///
unsigned int CLDPCDec::getK()
{
	return m_N-m_M;	
}

/// \brief	Returns the number of edges.
///
unsigned int CLDPCDec::getE()
{
	return m_E;	
}

/// \brief  Decoding of an LDPC code using defined method
/// \param	Lch			L-values received from the channel
/// \param	Lapp		a-posteriori L-values after decoding
/// \param	maxit		maximum number of iterations to perform
///
/// On termination the number of iterations performed is returned.
unsigned int CLDPCDec::decode(MSGTYPE *Lch, MSGTYPE *Lapp, unsigned int maxit)
{
	switch(m_method)
	{
		case DECODER_SPA:
			return decodeSPA(Lch, Lapp, maxit);
		case DECODER_MSA:
			return decodeMSA(Lch, Lapp, maxit);
		default:
			return decodeSPA(Lch, Lapp, maxit);
	}
}

/// \brief  Decoding of an LDPC code using the sum-product algorithm.
/// \param	Lch			L-values received from the channel
/// \param	Lapp		a-posteriori L-values after decoding
/// \param	maxit		maximum number of iterations to perform
///
/// This is the implementation of the sum-product algorithm for decoding
/// of an LDPC code. The decoder performs up to maxit iterations or
/// terminates if all check nodes are satisfied.
/// On termination the number of iterations performed is returned.
unsigned int CLDPCDec::decodeSPA(MSGTYPE *Lch, MSGTYPE *Lapp, unsigned int maxit)
{
	unsigned int	m, n, e, d;
	unsigned int	j, i, it;
	int						parity;
	MSGTYPE				x;
	MSGTYPE 	d1,d2,i1;
	bool					valid;
	
	if(m_clearmsg)
	{
		// initialize messages (in LLR domain)
		for(e=0;e<m_E;e++)
			m_c2v[e] = 0.0;
	}

	// perform iterations
	valid = false;
	for(it=0; it<maxit; it++)
	{
		// varnode processing
		d = 0;
		for(n=0; n<m_N; n++)
		{
			Lapp[n] = Lch[n];
			for(j=0; j<m_vardegree[n]; j++)
				Lapp[n] += m_c2v[d+j];

			for(j=0; j<m_vardegree[n]; j++)
			{
		    m_v2c[d+j]	= EXP(-CLIP(Lapp[n] - m_c2v[d+j], MAXLLR));
				m_hard[d+j]	= (Lapp[n]<=0);
			}
					
			d += m_vardegree[n];
		}
		
		// chknode processing
		valid 	= true;
		d				= 0;
		for(m=0; m<m_M; m++)
		{
			parity = 0;
			for(j=0; j<m_chkdegree[m]; j++)
				parity ^= m_hard[m_interleaver[d+j]];

			if(parity)
				valid = false;
			
			for(j=0; j<m_chkdegree[m]; j++)
			{
				x		= 0.0;
				
				for(i=0; i<m_chkdegree[m]; i++)
					if(j != i)
					{
					    i1 = m_v2c[m_interleaver[d+i]];
					    d1 = x+i1;
					    d2 = 1+x*i1;
					    x= d1/d2;
					}
				m_c2v[m_interleaver[d+j]] = -LOG(x);
			}
			d += m_chkdegree[m];
		}
    
 		if(valid)
			break;		
   
	}
	
	return it;	
}

/// \brief  Decoding of an LDPC code using the min-sum algorithm.
/// \param	Lch			L-values received from the channel
/// \param	Lapp		a-posteriori L-values after decoding
/// \param	maxit		maximum number of iterations to perform
///
/// This is the implementation of the min-sum algorithm for decoding
/// of an LDPC code. The decoder performs up to maxit iterations or
/// terminates if all check nodes are satisfied.
/// On termination the number of iterations performed is returned.
unsigned int CLDPCDec::decodeMSA(MSGTYPE *Lch, MSGTYPE *Lapp, unsigned int maxit)
{
	unsigned int	m, n, e, d, dc;
	unsigned int	j, i, it;
	int						parity, hd;
	bool					valid;
  
	if(m_clearmsg)
	{
		// initialize messages (in LLR domain)
		for(e=0;e<m_E;e++)
			m_c2v[e] = 0.0;
	}
	
	// perform iterations
	for(it=0; it<maxit; it++)
	{
		// varnode processing
		d = 0;
		for(n=0; n<m_N; n++)
		{
			Lapp[n] = Lch[n];
			for(j=0; j<m_vardegree[n]; j++)
				Lapp[n] += m_c2v[d+j];
			
      hd = HD(Lapp[n]);
			for(j=0; j<m_vardegree[n]; j++)
			{
        m_tmp[d+j]	= Lapp[n] - m_c2v[d+j];
				m_hard[d+j]	= hd;
			}
			
			d += m_vardegree[n];
		}
    
    for(e=0;e<m_E;e++)
      m_v2c[e] = CLIP(m_tmp[m_interleaver[e]], MAXLLR);
		
		// chknode processing
		valid 	= true;
		d				= 0;
		for(m=0; m<m_M; m++)
		{
      dc = m_chkdegree[m];
      
			parity = 0;
			for(j=0; j<dc; j++)
				parity ^= m_hard[m_interleaver[d+j]];
			
			if(parity)
				valid = false;
      
      // forward/backward algorithm for check node
      // init
      m_mf[0]    = FABS(m_v2c[d]);
      m_mb[dc-2] = FABS(m_v2c[d+dc-1]);
      m_sf[0]    = HD(m_v2c[d]);
      m_sb[dc-2] = HD(m_v2c[d+dc-1]);
      
      // forward
      for(j=1;j<(dc-1);j++)
      {
        m_mf[j]       = MIN(m_mf[j-1], FABS(m_v2c[d+j]));
        m_sf[j]       = m_sf[j-1] ^ HD(m_v2c[d+j]);
      }
      // backward
      for(i=dc-2;i>0;i--)
      {
        m_mb[i-1]  = MIN(m_mb[i], FABS(m_v2c[d+i]));
        m_sb[i-1]  = m_sb[i] ^ HD(m_v2c[d+i]);
      }
      // combine
      m_tmp[d]      = copysignf(m_mb[0], 1-2*(float)m_sb[0]);
      m_tmp[d+dc-1] = copysignf(m_mf[dc-2], 1-2*(float)m_sf[dc-2]);
      for(j=1;j<(dc-1);j++)
        m_tmp[d+j] = copysignf(MIN(m_mf[j-1], m_mb[j]), 1-2*(float)(m_sf[j-1]^m_sb[j]));
      
			d += dc;
		}
    
		if(valid)
			break;
    
    // post-processing
    if(it<m_corrveclen)
    {
      for(e=0;e<m_E;e++)
        m_c2v[e] *= m_corrvec[it];
    }
    
    for(e=0;e<m_E;e++)
      m_c2v[m_interleaver[e]] = m_tmp[e];
  }
	
	return it;	
}

/// \brief  Reads parity-check matrix from file.
/// \param	filename	name of the file
///
/// Reads a parity-check matrix from a file which is in our format.
int CLDPCDec::readDecoder(const char *filename)
{
	FILE          *file;
	unsigned int	n, m, e;

	file = fopen(filename, "r");
	if(file==NULL)
		return 0;
	
	fscanf(file, "%u %u %u", &m_N, &m_M, &m_E);
  
	allocateMem();
    
	for(n=0;n<m_N;n++)
		fscanf(file, "%u", &(m_vardegree[n]));
	for(m=0;m<m_M;m++)
		fscanf(file, "%u", &(m_chkdegree[m]));
	for(e=0;e<m_E;e++)
		fscanf(file, "%u", &(m_interleaver[e]));

	fclose(file);

  m_dcmax = 0;
  for(m=0;m<m_M;m++)
    m_dcmax = MAX(m_dcmax, m_chkdegree[m]);
  
  m_mf = (MSGTYPE*)malloc(m_dcmax*sizeof(MSGTYPE));
  m_mb = (MSGTYPE*)malloc(m_dcmax*sizeof(MSGTYPE));
  m_sf = (unsigned int*)malloc(m_dcmax*sizeof(unsigned int));
  m_sb = (unsigned int*)malloc(m_dcmax*sizeof(unsigned int));
  
	return 1;
}

/// \brief  Sets parameter clearmsg.
/// \param	clearmsg	new value
///
/// Sets parameter clearmsg.
void CLDPCDec::set_param_clearmsg(bool clearmsg)
{
	m_clearmsg = clearmsg;
}

/// \brief  Sets parameter method.
/// \param	method	new value
///
/// Sets parameter method.
void CLDPCDec::set_param_method(unsigned int method)
{
	m_method = method;
}

/// \brief  Sets correction vector for MSA post-processing
/// \param	len      number of correction vectors
/// \param  corrvec  pointer to vector
///
/// Sets parameter method.
void CLDPCDec::set_param_corrvec(unsigned int len, MSGTYPE *corrvec)
{
  unsigned int i;
  
  free(m_corrvec);
  m_corrvec = (MSGTYPE*)malloc(len * sizeof(MSGTYPE));
  m_corrveclen = len;
  
  for(i=0;i<m_corrveclen;i++)
    m_corrvec[i] = corrvec[i];
}

/// \brief	Determines whether the decoder is loaded
///
/// Determines whether the decoder is loaded
bool CLDPCDec::dec_loaded()
{
	return m_decloaded;
}

unsigned int CLDPCDec::encodeLDGM(unsigned int *info, unsigned int *parity)
{
	unsigned int	m, k, d;
	unsigned int	j;
  
  // varnode processing
  d = 0;
  for(k=0; k<getK(); k++)
  {
    for(j=0; j<m_vardegree[k]; j++)
      m_hard[d+j]	= info[k];
    
    d += m_vardegree[k];
  }
  
  // interleaving and chknode processing
  d				= 0;
  for(m=0; m<m_M; m++)
  {
    parity[m] = 0;
    for(j=0; j<(m_chkdegree[m]-1); j++)
      parity[m] ^= m_hard[m_interleaver[d+j]];
        
    d += m_chkdegree[m];
  }
  		
	return 1;	
}

unsigned int CLDPCDec::encodeRA(unsigned int *info, unsigned int *parity)
{
	unsigned int	m, k, dv, dc;
	unsigned int	j;
  
  // varnode processing
  dv = 0;
  for(k=0; k<getK(); k++)
  {
    for(j=0; j<m_vardegree[k]; j++)
      m_hard[dv+j]	= info[k];
    
    dv += m_vardegree[k];
  }
  
  // interleaving, chknode processing and accumulating
  dc				= 0;
  for(m=0; m<m_M; m++)
  {
    parity[m] = 0;
    for(j=0; j<(m_chkdegree[m]-1); j++)
      parity[m] ^= m_hard[m_interleaver[dc+j]];
    
    m_hard[dv+1] = parity[m];
    
    dv += 2;
    dc += m_chkdegree[m];
  }
  
	return 1;
}

