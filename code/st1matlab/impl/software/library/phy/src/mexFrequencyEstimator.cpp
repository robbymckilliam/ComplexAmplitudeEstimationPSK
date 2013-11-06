//
//  mexFrequencyEstimator.cpp
//  mexFrequencyEstimator
//
//  Created by Andre Pollok on 07/02/2013.
//  Copyright (c) 2013 ITR, University of South Australia. All rights reserved.
//

#include "mex.h"
#include "MatlabHelper.h"
#include "FrequencyEstimator.h"
#include <map>
#include <string>
#include <sstream> 

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

void Cleanup();

void DisplayHelp();

static std::map<std::string,FrequencyEstimator*> freqEstMap; // memory for estimators
static std::vector<complexT> rxSamp, txHatSamp; // memory for received and estimated transmit samples
static std::vector<VALTYPE> tSamp; // memory for sample times


/// \brief	Frees allocated memory.
///
void Cleanup()
{
  //delete all the estimators in the map.
  std::map<std::string,FrequencyEstimator*>::iterator itr;
  for(itr = freqEstMap.begin(); itr!=freqEstMap.end(); ++itr) delete itr->second;
  freqEstMap.clear();
  rxSamp.clear();
  txHatSamp.clear();
  tSamp.clear();
}

/// \brief  Displays help text
///
void DisplayHelp()
{
  mexPrintf("\nmexFrequencyEstimator MATLAB wrapper the coherent phase estimator\n");
  mexPrintf("\n(C) 2013 Andre Pollok, ITR/UniSA\n");
  mexPrintf("Compiled: %s %s\n\n", __DATE__, __TIME__);
  mexPrintf("Messages are represented as %s.\n\n", STRINGIZE_VALUE_OF(MSGTYPE));
}

/// \brief  Main function for MATLAB calls.
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
  mexAtExit(Cleanup);	
	
  // construct least-squares frequency estimator
  if (nrhs == 4){
    std::string name = mxArrayToString(prhs[0]);
    const unsigned int Nfft = (const unsigned int)mxGetScalar(prhs[1]);
    const VALTYPE Ts = (const VALTYPE)mxGetScalar(prhs[2]);
    const VALTYPE tauref = (const VALTYPE)mxGetScalar(prhs[3]);
    
    delete freqEstMap[name]; // delete item of same name if it exists
    freqEstMap[name] = new FrequencyEstimatorLeastSquares(Nfft,Ts,tauref);
    return;
  }
  
  else if (nrhs == 7) {
    std::string name = mxArrayToString(prhs[0]);
    const unsigned int Nfft = (const unsigned int)mxGetScalar(prhs[1]);
    const VALTYPE Ts = (const VALTYPE)mxGetScalar(prhs[2]);
    const VALTYPE tauref = (const VALTYPE)mxGetScalar(prhs[3]);
    const VALTYPE frmin = (const VALTYPE)mxGetScalar(prhs[4]);
    const VALTYPE frmax = (const VALTYPE)mxGetScalar(prhs[5]);
    const VALTYPE search_oversample = (const VALTYPE)mxGetScalar(prhs[6]);
    
    delete freqEstMap[name]; // delete item of same name if it exists
    freqEstMap[name] = new FrequencyEstimatorLeastSquaresWithRate(Nfft,Ts,tauref,frmin,frmax,search_oversample);
    return;
  }

  // run least-squares frequency estimator 
  else if ((nrhs == 2) | (nrhs == 3)){
    //get the estimator from the map
    std::string name = mxArrayToString(prhs[0]);
    // if the estimator hasn't been created, tell the user 
    if(freqEstMap.find(name) == freqEstMap.end()){
      mexPrintf("Estimator named '%s' has not been created \n\n", name.c_str());
      return;
    }
    //otherwise get the estimator
    FrequencyEstimator& est = *freqEstMap[name];
    
    //get the signal and run the estimator
    CopyMatlabComplexToVector(prhs[1], rxSamp);
    if (nrhs == 2) {
      est.estimate(rxSamp);
    } else if (nrhs == 3) {
      CopyMatlabComplexToVector(prhs[2], txHatSamp);
      est.estimate(rxSamp,txHatSamp);
    }
    est.refineEstimate();

    //mexPrintf("mex out chat = %f + %f i\n", real(est.complexGainEstimate()), imag(est.complexGainEstimate()));

    // set the estimator output
    if ( nlhs > 0 ) plhs[0] = mxCreateDoubleScalar(est.frequencyEstimate());
    if ( nlhs > 1 ) plhs[1] = ToMatlabComplex(est.complexGainEstimate());
    if ( nlhs > 2 ) plhs[2] = mxCreateDoubleScalar(est.frequencyRateEstimate());
    if ( nlhs > 3 ) plhs[3] = mxCreateDoubleScalar(est.objectiveFunctionValue());
    
    return;
  }

  // check basic parameter format
  else {
    DisplayHelp();
    return;
  }

}
