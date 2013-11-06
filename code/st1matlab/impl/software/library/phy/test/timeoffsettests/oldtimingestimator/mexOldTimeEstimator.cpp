//
//  mexTimingEstimator.cpp
//  TimingEstimator
//
//  Created by Andre Pollok on 31/05/2012.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#include "mex.h"
#include <iostream>
#include <complex>
#include "MatlabHelper.h"
#include "TimingEstimator.h"
#include <math.h>
#include <string>
#include <map>

using namespace std;

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

static std::map<std::string,PulseShape*> PulseShapeMap; // memory for pulse shapes
static std::map<std::string,TimingEstimator*> TimingEstimatorMap; // memory for timing estimators

void Cleanup();
void DisplayHelp();

void Cleanup()
{
  //delete all the pulse shapes in the map.
  std::map<std::string,PulseShape*>::iterator itr;
  for(itr = PulseShapeMap.begin(); itr!=PulseShapeMap.end(); ++itr) delete itr->second;
  PulseShapeMap.clear();

  //delete all the estimators in the map.
  std::map<std::string,TimingEstimator*>::iterator itr2;
  for(itr2 = TimingEstimatorMap.begin(); itr2!=TimingEstimatorMap.end(); ++itr2) delete itr2->second;
  TimingEstimatorMap.clear();
}

/// \brief  Displays help text
///
void DisplayHelp()
{
  mexPrintf("\nmexTimingEstimator   MATLAB wrapper for TimingEstimator class\n");
  mexPrintf("This mex-file enables the use of the TimingEstimator class from MATLAB\n");
  mexPrintf("\n(C) 2012 Andre Pollok, ITR/UniSA\n");
  mexPrintf("Compiled: %s %s\n\n", __DATE__, __TIME__);
  mexPrintf("Floating point numbers are represented as %s.\n\n", STRINGIZE_VALUE_OF(VALTYPE));  
  mexPrintf("Please refer to mexTimingEstimator_help.m for detailed help.\n\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mexAtExit(Cleanup);
  
  VALTYPE          taumin;
  VALTYPE          taumax;
  VALTYPE          *rxSampRe = NULL;
  VALTYPE          *rxSampIm = NULL;
  complexT         *rxSamp = NULL;
  unsigned int     numSamp;
  unsigned int     *D = NULL;
  unsigned int     numData;
  VALTYPE          *pilotRe = NULL;
  VALTYPE          *pilotIm = NULL;
  complexT         *pilot = NULL;
  unsigned int     *P = NULL;
  unsigned int     numPilots;
  
  unsigned int     i;
  VALTYPE          *tauhat = NULL;
  VALTYPE          *fval = NULL;
  mwSize           dims[2];
  TimingEstimator* tauest;
  
  // construct pulse shape
  if (nrhs == 5){
    std::string pulseshape = mxArrayToString(prhs[0]);
    VALTYPE rolloff, Tsymb, OS, Tpshape;
    rolloff = *((MSGTYPE*) mxGetData(prhs[1]));
    Tsymb = *((MSGTYPE*) mxGetData(prhs[2]));
    OS = *((MSGTYPE*) mxGetData(prhs[3]));
    Tpshape = *((MSGTYPE*) mxGetData(prhs[4]));
    delete PulseShapeMap[pulseshape];
    PulseShapeMap[pulseshape] = new RootRaisedCosine(rolloff, Tsymb, OS, Tpshape);
    
  // construct time offset estimator
  } else if (nrhs == 7 || nrhs == 8 || nrhs == 9) {

    // assign input parameters
    std::string estimatorName = mxArrayToString(prhs[0]);
    std::string pulseshapeName = mxArrayToString(prhs[1]);
    // if the pulse shape hasn't been created, tell the user 
    if(PulseShapeMap.find(pulseshapeName) == PulseShapeMap.end()){
      mexPrintf("Pulse shape named '%s' has not been created \n\n", pulseshapeName.c_str());
      return;
    }
    
    numData = (unsigned int)mxGetN(prhs[2]);
    D = new unsigned int [numData];
    GetIntArray(prhs[2], D);
    
    numPilots = (unsigned int)mxGetN(prhs[3]);
    P = new unsigned int [numPilots];
    GetIntArray(prhs[3], P);
    
    pilotRe = new VALTYPE [numPilots];
    pilotIm = new VALTYPE [numPilots];
    pilot   = new complexT [numPilots];
    GetFloatArray(prhs[4], pilotRe, false);
    GetFloatArray(prhs[4], pilotIm, true);
    for (i=0; i<numPilots; i++) {
      pilot[i] = complexT (pilotRe[i],pilotIm[i]);
    }
    
    taumin = (VALTYPE)mxGetScalar(prhs[5]);
    taumax = (VALTYPE)mxGetScalar(prhs[6]);
    
    delete TimingEstimatorMap[estimatorName]; //delete element with this name if it exists
    if (nrhs <= 8) {
      TimingEstimatorMap[estimatorName] = new TimingEstimator(taumin,taumax,PulseShapeMap[pulseshapeName],D,numData,pilot,P,numPilots);
    } else if (nrhs == 9) {
      unsigned int a;
      a = (unsigned int)mxGetScalar(prhs[8]);
      TimingEstimatorMap[estimatorName] = new TimingEstimator(taumin,taumax,PulseShapeMap[pulseshapeName],D,numData,pilot,P,numPilots,a);
    }
    if (nrhs >= 8) {
      // set estimator tolerance
      VALTYPE tol;
      tol = (VALTYPE)mxGetScalar(prhs[7]);
      TimingEstimatorMap[estimatorName]->setTolerance(tol);
    }

    // clean up
    delete[] pilot;
    delete[] pilotIm;
    delete[] pilotRe;
    delete[] P;
    delete[] D;
    
  // run estimator
  } else if (nrhs == 2) {
    
    // assign input parameters
    std::string estimatorName = mxArrayToString(prhs[0]);
    // if the pulse shape hasn't been created, tell the user 
    if(TimingEstimatorMap.find(estimatorName) == TimingEstimatorMap.end()){
      mexPrintf("Timing estimator named '%s' has not been created \n\n", estimatorName.c_str());
      return;
    } else {
      //otherwise get the estimator
      tauest = TimingEstimatorMap[estimatorName];
    }
    
    numSamp  = (unsigned int)mxGetN(prhs[1]); 
    rxSampRe = new VALTYPE [numSamp];
    rxSampIm = new VALTYPE [numSamp];
    rxSamp   = new complexT [numSamp];
    GetFloatArray(prhs[1], rxSampIm, true);
    GetFloatArray(prhs[1], rxSampRe, false);
    for (i=0; i<numSamp; i++) {
      rxSamp[i] = complexT (rxSampRe[i],rxSampIm[i]);
    }
    
    // Allocate memory for output parameters
    dims[0] = 1;
    dims[1] = 1;
    plhs[0]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[1]), mxREAL);
    plhs[1]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[1]), mxREAL);
    // get pointer to outputs
    tauhat = (VALTYPE*)mxGetPr(plhs[0]);
    fval   = (VALTYPE*)mxGetPr(plhs[1]);
    
    // run estimator
    tauest->runEstimator(rxSamp,numSamp);
    *tauhat = tauest->getTimingEstimate();
    *fval   = tauest->getFuncValue();
  
// DEBUGGING  
/*    // Yk
    unsigned int c = 4;
    unsigned int K = floor((taumax-taumin)/(T/OS/c))+1;
    dims[1] = K;
    complexT* Yk = tauest->getYk();
    VALTYPE* YkRe;
    VALTYPE* YkIm;
    YkRe = new VALTYPE [K];
    YkIm = new VALTYPE [K];
    for (i=0; i<K; i++) {
      YkRe[i] = real(Yk[i]);
      YkIm[i] = imag(Yk[i]);
    }
    plhs[2]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[0]), mxCOMPLEX);
    SetFloatArray(YkRe, plhs[2], false);
    SetFloatArray(YkIm, plhs[2], true);
    delete[] YkRe;
    delete[] YkIm;
    // Zk
    VALTYPE* Zk = tauest->getZk();
    plhs[3]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[0]), mxREAL);
    SetFloatArray(Zk, plhs[3], false);
    // bk
//    K = 1905;
//    K = 67664;
//    K = 1649;
    K = 321;
    dims[1] = K;
    complexT* bk = tauest->getbkpub();
    VALTYPE* bkRe;
    VALTYPE* bkIm;
    bkRe = new VALTYPE [K];
    bkIm = new VALTYPE [K];
    for (i=0; i<K; i++) {
      bkRe[i] = real(bk[i]);
      bkIm[i] = imag(bk[i]);
    }
    plhs[4]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[0]), mxCOMPLEX);
    SetFloatArray(bkRe, plhs[4], false);
    SetFloatArray(bkIm, plhs[4], true);
    delete[] bkRe;
    delete[] bkIm;
/*  
    // gk
    K = 161;
    dims[1] = K;
    VALTYPE* gk = tauest->getgk();
    plhs[5]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[0]), mxREAL);
    SetFloatArray(gk, plhs[5], false);
    // zk
    K = c*numSamp;
    dims[1] = K;
    complexT* zk = tauest->getzk();
    VALTYPE* zkRe;
    VALTYPE* zkIm;
    zkRe = new VALTYPE [K];
    zkIm = new VALTYPE [K];
    for (i=0; i<K; i++) {
      zkRe[i] = real(zk[i]);
      zkIm[i] = imag(zk[i]);
    }
    plhs[6]	= mxCreateNumericArray(2, dims, mxGetClassID(prhs[0]), mxCOMPLEX);
    SetFloatArray(zkRe, plhs[6], false);
    SetFloatArray(zkIm, plhs[6], true);
    delete[] zkRe;
    delete[] zkIm;
*/
// END DEBUGGING  
  
    SetFloatArray(tauhat, plhs[0], false);
    SetFloatArray(fval,   plhs[1], false);
    
    // clean up
    delete[] rxSamp;
    delete[] rxSampIm;
    delete[] rxSampRe;
    
    return;
    
  } else {
    DisplayHelp();
  }
}

