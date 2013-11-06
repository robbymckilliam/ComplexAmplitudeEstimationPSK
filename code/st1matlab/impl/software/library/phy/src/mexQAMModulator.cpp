//
//  mexQAMModulator.cpp
//  QAMModulator
//
//  Created by Gottfried Lechner on 16/03/12.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#include "mex.h"
#include "QAMModulator.h"
#include <string.h>
#include <math.h>
#include "MatlabHelper.h"

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

void Cleanup();
void DisplayHelp();

/// \brief	defines the maximum number of modulator instances
///
#define	MAXMOD	100

static CQAMModulator		**pMod= NULL;

/// \brief	Frees allocated memory.
///
void Cleanup()
{
	unsigned int	i;
  
	for(i=0;i<MAXMOD;i++)
		delete(pMod[i]);
	free(pMod);
	pMod = NULL;
}

/// \brief  Displays help text
///
void DisplayHelp()
{
  mexPrintf("\nmexQAMModulator   MATLAB wrapper for QAMModulator class\n");
  mexPrintf("This mex-file enables the use of the QAMModulator class from MATLAB\n");
  mexPrintf("\n(C) 2012 Gottfried Lechner, ITR/UniSA\n");
  mexPrintf("Compiled: %s %s\n\n", __DATE__, __TIME__);
  mexPrintf("Messages are represented as %s.\n\n", STRINGIZE_VALUE_OF(MSGTYPE));  
  mexPrintf("Use mexQAMModulator(id, command)\n\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char					command[100], parameter[100];
	unsigned int	i, ID, n;
  mwSize        dims[2];
  
	mexAtExit(Cleanup);	
  
	if(!pMod)
	{
    pMod = (CQAMModulator**)malloc(MAXMOD * sizeof(CQAMModulator*));
		for(i=0;i<MAXMOD;i++)
			pMod[i] = NULL;
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
  if(ID>=MAXMOD)
    mexErrMsgTxt("ID exceeds MAXMOD.");
  
	if(!pMod[ID])
		pMod[ID] = new(CQAMModulator);
  
	mxGetString(prhs[1], command, 100);
  
	if(strcmp(command, "init")==0)
	{
    unsigned int  numsections, numsymbols;
    MSGTYPE       *SymReal, *SymImag;
    
    // get parameters and pointers to trellis structure
    numsections = (unsigned int)mxGetScalar(prhs[2]);
    numsymbols  = (unsigned int)mxGetN(prhs[3]);
    
    SymReal     = (MSGTYPE*)malloc(numsymbols * sizeof(MSGTYPE));
    SymImag     = (MSGTYPE*)malloc(numsymbols * sizeof(MSGTYPE));
    
    GetFloatArray(prhs[3], SymReal, false);
    GetFloatArray(prhs[3], SymImag, true);    
        
    pMod[ID]->init(numsections, numsymbols, SymReal, SymImag);
    
    free(SymReal);
    free(SymImag);
	}
	else if(strcmp(command, "set-ch")==0)
	{
    MSGTYPE  *yReal, *yImag, *sigma2;
    
    // check if initialised
    if(pMod[ID]->getlength()==0)
      mexErrMsgTxt("Modulator not initialised");
    
    // check length
    if(mxGetN(prhs[2]) != pMod[ID]->getlength())
      mexErrMsgTxt("set-ch: wrong length");

    yReal     = (MSGTYPE*)malloc(mxGetN(prhs[2]) * sizeof(MSGTYPE));
    yImag     = (MSGTYPE*)malloc(mxGetN(prhs[2]) * sizeof(MSGTYPE));
    sigma2    = (MSGTYPE*)malloc(mxGetN(prhs[2]) * sizeof(MSGTYPE));
    
    GetFloatArray(prhs[2], yReal, false);
    GetFloatArray(prhs[2], yImag, true);    
    
    if(mxGetN(prhs[3])==1)
    {
      for(n=0;n<mxGetN(prhs[2]);n++)
        sigma2[n] = (MSGTYPE)mxGetScalar(prhs[3]);
    }
    else
    {
      GetFloatArray(prhs[3], sigma2, false);
    }

    pMod[ID]->set_ch(yReal, yImag, sigma2);
    
    free(yReal);
    free(yImag);
    free(sigma2);
	}
	else if(strcmp(command, "set-apri")==0)
	{
    MSGTYPE  *La;

    // check if initialised
    if(pMod[ID]->getlength()==0)
      mexErrMsgTxt("Modulator not initialised");    
    
    // check length
    if(mxGetN(prhs[2]) != pMod[ID]->getlength())
      mexErrMsgTxt("set-apri: wrong length");

    La = (MSGTYPE*)malloc(pMod[ID]->getlength()*pMod[ID]->getbitspersymbol() * sizeof(MSGTYPE));
    
    GetFloatArray(prhs[2], La, false);

    pMod[ID]->set_apri(La);
    
    free(La);
	}
	else if(strcmp(command, "modulate")==0)
	{
    MSGTYPE       *yReal, *yImag;
    unsigned int  *data;

    // check if initialised
    if(pMod[ID]->getlength()==0)
      mexErrMsgTxt("Modulator not initialised");
    
    // check length
    if(mxGetN(prhs[2]) != pMod[ID]->getlength()*pMod[ID]->getbitspersymbol())
      mexErrMsgTxt("modulate: wrong length");
    
    data    = (unsigned int*)malloc(mxGetN(prhs[2]) * sizeof(unsigned int));
    yReal   = (MSGTYPE*)malloc(mxGetN(prhs[2]) * sizeof(MSGTYPE));
    yImag   = (MSGTYPE*)malloc(mxGetN(prhs[2]) * sizeof(MSGTYPE));
    
    GetIntArray(prhs[2], data);
    
    dims[0] = 1;
    dims[1] = pMod[ID]->getlength();
    plhs[0]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[2]), mxCOMPLEX);

		pMod[ID]->modulate(data, yReal, yImag);
    
    SetFloatArray(yReal, plhs[0], false);
    SetFloatArray(yImag, plhs[0], true);

    free(data);
    free(yReal);
    free(yImag);
 	}
	else if(strcmp(command, "demodulate")==0)
	{
    MSGTYPE   *Lout;

    // check if initialised
    if(pMod[ID]->getlength()==0)
      mexErrMsgTxt("Modulator not initialised");
    
    dims[0] = pMod[ID]->getbitspersymbol();
    dims[1] = pMod[ID]->getlength();
    
    plhs[0]	= mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    
    Lout    = (MSGTYPE*)malloc(pMod[ID]->getlength()*pMod[ID]->getbitspersymbol() * sizeof(MSGTYPE));
    
		pMod[ID]->demodulate(Lout);
    
    SetFloatArray(Lout, plhs[0]);
    
    free(Lout);
 	}
	else if(strcmp(command, "remodulate")==0)
	{
    MSGTYPE  *yReal, *yImag, *ResVar;

    // check if initialised
    if(pMod[ID]->getlength()==0)
      mexErrMsgTxt("Modulator not initialised");
    
    dims[0] = 1;
    dims[1] = pMod[ID]->getlength();
    plhs[0]	= mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);    
    
    yReal   = (MSGTYPE*)malloc(pMod[ID]->getlength() * sizeof(MSGTYPE));
    yImag   = (MSGTYPE*)malloc(pMod[ID]->getlength() * sizeof(MSGTYPE));
    
    if(nlhs>1)
    {
      dims[0] = 1;
      dims[1] = pMod[ID]->getlength();      
      plhs[1]	= mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);    
      ResVar  = (MSGTYPE*)malloc(pMod[ID]->getlength() * sizeof(MSGTYPE)); 
      pMod[ID]->remodulate(yReal, yImag, ResVar);
      SetFloatArray(ResVar, plhs[1], false);
    }
    else
    {
      ResVar = NULL;
      pMod[ID]->remodulate(yReal, yImag, NULL);
    }
    
    SetFloatArray(yReal, plhs[0], false);
    SetFloatArray(yImag, plhs[0], true);    
    
    free(yReal);
    free(yImag);
    free(ResVar);
  }
	else if(strcmp(command, "configure")==0)
	{
		mxGetString(prhs[2], parameter, 100);
  }
}

