//
//  mexQAMModulatorV2.cpp
//  QAMModulator
//
//  Created by Gottfried Lechner on 09/02/13.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#include "mex.h"
#include "MatlabHelper.h"
#include "QAMModulatorV2.h"
#include <string.h>
#include <iostream>
#include <math.h>
#include <complex>


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
  mexPrintf("\nmexQAMModulator   MATLAB wrapper for QAMModulator class Version 2\n");
  mexPrintf("This mex-file enables the use of the QAMModulator class from MATLAB\n");
  mexPrintf("\n(C) 2013 Gottfried Lechner, ITR/UniSA\n");
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
    
    // get parameters and pointers to trellis structure
    numsections = (unsigned int)mxGetScalar(prhs[2]);
    numsymbols  = (unsigned int)mxGetN(prhs[3]);
        
    pMod[ID]->init(numsections, numsymbols, GetComplexArray(prhs[3]));
	}
	else if(strcmp(command, "set-ch")==0)
	{
    // check if initialised
    if(pMod[ID]->getlength()==0)
      mexErrMsgTxt("Modulator not initialised");
    
    // check length
    if(mxGetN(prhs[2]) != pMod[ID]->getlength())
      mexErrMsgTxt("set-ch: wrong length");
    
    if(mxGetN(prhs[3])==1)
      pMod[ID]->set_ch(GetComplexArray(prhs[2]), (MSGTYPE)mxGetScalar(prhs[3]));
    else
      pMod[ID]->set_ch(GetComplexArray(prhs[2]), GetFloatArray(prhs[3]));
	}
	else if(strcmp(command, "set-apri")==0)
	{
    // check if initialised
    if(pMod[ID]->getlength()==0)
      mexErrMsgTxt("Modulator not initialised");    
    
    // check length
    if(mxGetN(prhs[2]) != pMod[ID]->getlength())
      mexErrMsgTxt("set-apri: wrong length");

    pMod[ID]->set_apri(GetFloatArray(prhs[2]));
 	}
	else if(strcmp(command, "modulate")==0)
	{
    unsigned int  *data;

    // check if initialised
    if(pMod[ID]->getlength()==0)
      mexErrMsgTxt("Modulator not initialised");
    
    // check length
    if(mxGetN(prhs[2]) != pMod[ID]->getlength()*pMod[ID]->getbitspersymbol())
      mexErrMsgTxt("modulate: wrong length");
    
    data    = (unsigned int*)malloc(mxGetN(prhs[2]) * sizeof(unsigned int));
    
    GetIntArray(prhs[2], data);    
    std::vector<complex>  y(pMod[ID]->getlength());
    
		pMod[ID]->modulate(data, y);
    
    plhs[0] = ToMatlabComplex(y);

    free(data);
 	}
	else if(strcmp(command, "demodulate")==0)
	{
    // check if initialised
    if(pMod[ID]->getlength()==0)
      mexErrMsgTxt("Modulator not initialised");
    
    dims[0] = pMod[ID]->getbitspersymbol();
    dims[1] = pMod[ID]->getlength();
    
    plhs[0]	= mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
   
    std::vector<MSGTYPE>  Lout(pMod[ID]->getlength()*pMod[ID]->getbitspersymbol());
            
		pMod[ID]->demodulate(Lout);
    
    SetFloatArray(&Lout[0], plhs[0]);
 	}
	else if(strcmp(command, "remodulate")==0)
	{
    // check if initialised
    if(pMod[ID]->getlength()==0)
      mexErrMsgTxt("Modulator not initialised");
    
    std::vector<complex>  y(pMod[ID]->getlength());
    std::vector<MSGTYPE>  ResVar(pMod[ID]->getlength());
    
    pMod[ID]->remodulate(y, ResVar);
    
    plhs[0] = ToMatlabComplex(y);

    if(nlhs>1)
    {
      dims[0] = 1;
      dims[1] = pMod[ID]->getlength();
      plhs[1]	= mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      SetFloatArray(&ResVar[0], plhs[1]);
    }
  }
	else if(strcmp(command, "configure")==0)
	{
		mxGetString(prhs[2], parameter, 100);
  }
}

