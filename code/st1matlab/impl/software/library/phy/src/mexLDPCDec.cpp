#include "mex.h"
#include "LDPCDec.h"
#include <string.h>
#include "MatlabHelper.h"

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

void Cleanup();

void DisplayHelp();

/// \brief	defines the maximum number of LDPC instances
///
#define	MAXLDPC	200

static CLDPCDec		**pLDPC	= NULL;

/// \brief	Frees allocated memory.
///
void Cleanup()
{
	unsigned int	i;

	for(i=0;i<MAXLDPC;i++)
		delete(pLDPC[i]);
	free(pLDPC);
	pLDPC = NULL;
}

/// \brief  Displays help text
///
void DisplayHelp()
{
  mexPrintf("\nmexLDPCDec   MATLAB wrapper for LDPC class\n");
  mexPrintf("This mex-file enables the use of the LDPC class from MATLAB\n");
  mexPrintf("\n(C) 2009-2012 Gottfried Lechner, ITR/UniSA\n");
  mexPrintf("Compiled: %s %s\n\n", __DATE__, __TIME__);
  mexPrintf("Messages are represented as %s.\n\n", STRINGIZE_VALUE_OF(MSGTYPE));
  mexPrintf("Use mexLDPCDec(id, command)\n\n");
}

/// \brief  Main function for MATLAB calls.
/// \param	nlhs		number of left-hand-side arguments
/// \param	plhs		array of left-hand-side arguments
/// \param	nrhs		number of right-hand-side arguments
/// \param	prhs		array of right-hand-side arguments
///
///	The first argument is the ID of the LDPC instance (starting from 0).
/// This allows to use many LDPC codes in parallel. The maximum number of
/// parallel codes is defined as MAXLDPC (default=100).
///
/// Depending on the second argument, different tasks are performed and the
/// remaining arguments have a different meaning.
/// Possible commands (second argument) and further arguments:
/// - mexLDPCDec(ID, 'readdecoder', filename)
/// - (Lapp, it, Isyn)	= mexLDPCDec(ID, 'decode', Lch, maxit)
/// - (N, M)            = mexLDPCDec(ID, 'dimensions')
/// - mexLDPCDec(ID, 'configure', parameter, value)
///   - 'clearmsg'\n
///		If true, the decoder clears the messages from check to variable nodes.\n
///		If false, the decoder uses the messages from the previous call.
///   - 'method'\n
///   Defines the decoding method: 0 ... sum-product (default)
///                                1 ... min-sum
///   - 'corrvec'\n
///   Sets the correction vector for MSA decoding with post-processing
/// - parity = mexLDPCDec(ID, 'encodeLDGM', info)
///   Assumes that the code represents an LDGM code and computes parity bits.
/// - parity = mexLDPCDec(ID, 'encodeRA', info)
///   Assumes that the code represents an RA code and computes parity bits.
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char					command[100];
	char					parameter[100];
	char					filename[255];
	double				*ptr;
	unsigned int	it, maxit;
	unsigned int	i, ID;
  mwSize        dims[2];
  
	mexAtExit(Cleanup);	

	if(!pLDPC)
	{
		pLDPC = (CLDPCDec**)malloc(MAXLDPC * sizeof(CLDPCDec*));
		for(i=0;i<MAXLDPC;i++)
			pLDPC[i] = NULL;
	}

	// check basic parameter format
	if(nrhs<2)
  {
    DisplayHelp();
    return;
  }
	if(mxIsEmpty(prhs[0]) || mxIsEmpty(prhs[1]))
  {
    DisplayHelp();
    return;
  }
	   
	ID = (unsigned int)mxGetScalar(prhs[0]);
  if(ID>=MAXLDPC)
    mexErrMsgTxt("ID exceeds MAXLDPC.");

	if(!pLDPC[ID])
		pLDPC[ID] = new(CLDPCDec);

	mxGetString(prhs[1], command, 100);

	if(strcmp(command, "readdecoder")==0)
	{
		mxGetString(prhs[2], filename, 255);
		if(!(pLDPC[ID]->readDecoder(filename)))
			 mexErrMsgTxt("Unable to read decoder definition.");
	}
	else if(strcmp(command, "decode")==0)
	{
    MSGTYPE   *Lch, *Lapp;
    
    if(!pLDPC[ID]->dec_loaded())
      mexErrMsgTxt("Need to load decoder first.");
    
    dims[0] = 1;
    dims[1] = pLDPC[ID]->getN();
    plhs[0]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[2]), mxREAL);
    
    Lch     = (MSGTYPE*)malloc(pLDPC[ID]->getN() * sizeof(MSGTYPE));
    Lapp		= (MSGTYPE*)malloc(pLDPC[ID]->getN() * sizeof(MSGTYPE));
    
    GetFloatArray(prhs[2], Lch);
    
		maxit		= (unsigned int)mxGetScalar(prhs[3]);
    
		it			= pLDPC[ID]->decode(Lch, Lapp, maxit);
    
    SetFloatArray(Lapp, plhs[0]);
    free(Lapp);
    free(Lch);
    
		if(nlhs==1)
			return;
		
    dims[0] = 1;
    dims[1] = 1;
    plhs[1]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[2]), mxREAL);
    
    SetIntArray(&it, plhs[1]);
    
    if(nlhs==2)
      return;
    
    dims[0] = 1;
    dims[1] = 1;
    plhs[2]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[2]), mxREAL);
    
    MSGTYPE Isyn;
    Isyn    = pLDPC[ID]->syndromeInformation();
    
    SetFloatArray(&Isyn, plhs[2]);
	}
	else if(strcmp(command, "dimensions")==0)
	{
		plhs[0]	= mxCreateDoubleMatrix(1,1,mxREAL);
		plhs[1]	= mxCreateDoubleMatrix(1,1,mxREAL);

		ptr			= mxGetPr(plhs[0]);
		ptr[0]	= pLDPC[ID]->getN();

		ptr			= mxGetPr(plhs[1]);
		ptr[0]	= pLDPC[ID]->getM();
	}
	else if(strcmp(command, "configure")==0)
	{
		mxGetString(prhs[2], parameter, 100);

		if(strcmp(parameter, "clearmsg")==0)
			pLDPC[ID]->set_param_clearmsg((bool)mxGetScalar(prhs[3]));
		else if(strcmp(parameter, "method")==0)
			pLDPC[ID]->set_param_method((unsigned int)mxGetScalar(prhs[3]));
    else if(strcmp(parameter, "corrvec")==0)
    {
      MSGTYPE *corrvec;
      corrvec = (MSGTYPE*)malloc((unsigned int)mxGetN(prhs[3]) * sizeof(MSGTYPE));
      
      GetFloatArray(prhs[3], corrvec);
      
      pLDPC[ID]->set_param_corrvec((unsigned int)mxGetN(prhs[3]), corrvec);
      free(corrvec);
    }
	}
	else if(strcmp(command, "decloaded")==0)
	{
		plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
		ptr			= mxGetPr(plhs[0]);
		ptr[0]	= (double)(pLDPC[ID]->dec_loaded());
	}
  else if(strcmp(command, "encodeLDGM")==0)
	{
    if(!pLDPC[ID]->dec_loaded())
      mexErrMsgTxt("Need to load decoder first.");    
    
    unsigned int *info, *parity;
    
    info		= (unsigned int*)malloc(pLDPC[ID]->getN() * sizeof(unsigned int));
		parity	= (unsigned int*)malloc(pLDPC[ID]->getM() * sizeof(unsigned int));

    dims[0] = 1;
    dims[1] = pLDPC[ID]->getM();
    plhs[0]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[2]), mxREAL);
        
    GetIntArray(prhs[2], info);
        
    pLDPC[ID]->encodeLDGM(info, parity);
    
    SetIntArray(parity, plhs[0]);
    
		free(info);
		free(parity);
	}
  else if(strcmp(command, "encodeRA")==0)
	{
    if(!pLDPC[ID]->dec_loaded())
      mexErrMsgTxt("Need to load decoder first.");    
    
    unsigned int *info, *parity;
    
		info		= (unsigned int*)malloc(pLDPC[ID]->getK() * sizeof(unsigned int));
		parity	= (unsigned int*)malloc(pLDPC[ID]->getM() * sizeof(unsigned int));
    
    dims[0] = 1;
    dims[1] = pLDPC[ID]->getM();
    plhs[0]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[2]), mxREAL);    
    
    GetIntArray(prhs[2], info);
        
		pLDPC[ID]->encodeRA(info, parity);
    
    SetIntArray(parity, plhs[0]);
    
		free(info);
		free(parity);
	}
}
