//
//  mexTimingEstimator.cpp
//  TimingEstimator
//
//  Created by Andre Pollok on 31/05/2012.
//  Modifed by Robby McKilliam 16/2/2012.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#include "mex.h"
#include <iostream>
#include <complex>
#include "MatlabHelper.h"
#include "FinitePulse.h"
#include "TimeOffsetEstimator.h"
#include "FastAlgorithm.h"
#include "TabulatedFastAlgorithm.h"
#include "OerderMeyerMassey.h"
#include <math.h>
#include <string>
#include <map>

#ifndef VALTYPE
#define VALTYPE double
#endif
#ifndef MSGTYPE
#define MSGTYPE VALTYPE
#endif

using namespace std;
using namespace TimeOffset;

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

static std::map<std::string,TimeOffsetEstimator*> TimeOffsetEstimatorMap; // memory for timing estimators

void Cleanup();
void DisplayHelp();

void Cleanup()
{
  //delete all the estimators in the map.
  std::map<std::string,TimeOffsetEstimator*>::iterator itr2;
  for(itr2 = TimeOffsetEstimatorMap.begin(); itr2!=TimeOffsetEstimatorMap.end(); ++itr2) delete itr2->second;
  TimeOffsetEstimatorMap.clear();
}

/// \brief  Displays help text
///
void DisplayHelp()
{
  mexPrintf("\nmexTimeOffsetEstimator   MATLAB wrapper for TimeOffsetEstimator class\n");
  mexPrintf("This mex-file enables the use of the TimeOffsetEstimator class from MATLAB\n");
  mexPrintf("\n(C) 2012 Andre Pollok, ITR/UniSA\n");
  mexPrintf("Compiled: %s %s\n\n", __DATE__, __TIME__);
  mexPrintf("Floating point numbers are represented as %s.\n\n", STRINGIZE_VALUE_OF(VALTYPE));  
  mexPrintf("Please refer to mexTimeOffsetEstimator_help.m for detailed help.\n\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mexAtExit(Cleanup);

  if (nrhs == 12 || nrhs == 13) { //construct st1 recursive estimator 

    // assign input parameters
    std::string estimatorName = mxArrayToString(prhs[0]);
    std::vector<int> P = GetIndexArray(prhs[1]); //mexPrintf("P = "); for(unsigned int i = 0; i < P.size(); i++) mexPrintf("%i, ", P[i]);
    std::vector<int> D = GetIndexArray(prhs[2]); //mexPrintf("D = "); for(unsigned int i = 0; i < D.size(); i++) mexPrintf("%i, ", D[i]);
    std::vector<complex<VALTYPE> > pilots = GetComplexArray(prhs[3]);    
    VALTYPE T = (VALTYPE) mxGetScalar(prhs[4]); //symbol period
    VALTYPE Ts = (VALTYPE) mxGetScalar(prhs[5]); //sample period
    unsigned int p = (unsigned int) round( (VALTYPE) mxGetScalar(prhs[6]) );
    unsigned int q = (unsigned int) round( (VALTYPE) mxGetScalar(prhs[7]) );
    VALTYPE taumin = (VALTYPE) mxGetScalar(prhs[8]);
    VALTYPE taumax = (VALTYPE) mxGetScalar(prhs[9]);
    VALTYPE rolloff = (VALTYPE) mxGetScalar(prhs[10]);
    VALTYPE duration = (VALTYPE) mxGetScalar(prhs[11]);
    VALTYPE brenttol = (nrhs == 13 ? (VALTYPE)mxGetScalar(prhs[12]) : 1e-6); //default tolerance is 1e-6
    int tablesize = 50000; //tablesize is 50000 per symbol
    unsigned int c = 16; //c is always 16 (12 is about the smallest, this leaves some room).
    TruncatedRootRaisedCosine rrcpulse(T, rolloff, duration, true);
    delete TimeOffsetEstimatorMap[estimatorName]; //delete item of same name if it exists 
    //TimeOffsetEstimatorMap[estimatorName] = new Recursive<TruncatedRootRaisedCosine>(P,D,pilots,rrcpulse,T,Ts,taumin,taumax,c,p,q,brenttol);
    TimeOffsetEstimatorMap[estimatorName] = new TabulatedRecursive<TruncatedRootRaisedCosine>(P,D,pilots,rrcpulse,T,Ts,taumin,taumax,c,p,q,tablesize,brenttol);

  } else if (nrhs == 10) {  //construct st2 OerderMeyer estimator 
    
    // assign input parameters
    std::string estimatorName = mxArrayToString(prhs[0]);
    std::vector<int> P = GetIndexArray(prhs[1]); //mexPrintf("P = "); for(unsigned int i = 0; i < P.size(); i++) mexPrintf("%i, ", P[i]);
    std::vector<int> D = GetIndexArray(prhs[2]); //mexPrintf("D = "); for(unsigned int i = 0; i < D.size(); i++) mexPrintf("%i, ", D[i]);
    std::vector<complex<VALTYPE> > pilots = GetComplexArray(prhs[3]);    
    VALTYPE T = (VALTYPE) mxGetScalar(prhs[4]); //symbol period
    VALTYPE Ts = (VALTYPE) mxGetScalar(prhs[5]); //sample period
    VALTYPE taumin = (VALTYPE) mxGetScalar(prhs[6]);
    VALTYPE taumax = (VALTYPE) mxGetScalar(prhs[7]);
    VALTYPE rolloff = (VALTYPE) mxGetScalar(prhs[8]);
    VALTYPE duration = (VALTYPE) mxGetScalar(prhs[9]);
    int tablesize = 4000; //tablesize is 4000 per symbol
    TruncatedRootRaisedCosine rrcpulse(T, rolloff, duration, true);
    delete TimeOffsetEstimatorMap[estimatorName]; //delete item of same name if it exists 
    TimeOffsetEstimatorMap[estimatorName] = new TabulatedOerderMeyerAndMassey<TruncatedRootRaisedCosine>(P,D,pilots,rrcpulse,T,Ts,taumin,taumax,tablesize);
    
  // run estimator
  } else if (nrhs == 2) {
    
    TimeOffsetEstimator* tauest;
    // assign input parameters
    std::string estimatorName = mxArrayToString(prhs[0]);
    // if the pulse shape hasn't been created, tell the user 
    if(TimeOffsetEstimatorMap.find(estimatorName) == TimeOffsetEstimatorMap.end()){
      mexPrintf("Timing estimator named '%s' has not been created \n\n", estimatorName.c_str());
      return;
    } else { //otherwise get the estimator
      tauest = TimeOffsetEstimatorMap[estimatorName];
    }
    //get recieved signal
    std::vector<complex<VALTYPE> > rx = GetComplexArray(prhs[1]);    
    
    // run estimator
    VALTYPE tauhat = tauest->estimate(rx);
    
    //set the estimator output 
    if( nlhs > 0 ) plhs[0] = mxCreateDoubleScalar(tauhat);
    //if( nlhs > 1 ) plhs[1] = mxCreateDoubleScalar(tauest->SS(tauhat)); //maximum of objective function if second output argument requrested.
    
  } else {
    DisplayHelp();
  }

}

