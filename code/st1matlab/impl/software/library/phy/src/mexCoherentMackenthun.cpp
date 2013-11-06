#include "mex.h"
#include "CoherentMackenthun.h"
#include "MatlabHelper.h"
#include <map>
#include <string>
#include <sstream> 

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

void Cleanup();

void DisplayHelp();

static std::map<std::string,PhaseEstimator*> cmap; //memory for estimators
static std::vector<complex> y; //memory for input data
static std::vector<complex> p; //memory for pilots

#ifdef _WIN32
#define log2(x) 			(log(x)/log((MSGTYPE)2.0))
#define round(x)			(floor((x)+0.5))
#endif


/// \brief	Frees allocated memory.
///
void Cleanup()
{
  //delete all the estimators in the map.
  std::map<std::string,PhaseEstimator*>::iterator itr;
  for(itr = cmap.begin(); itr!=cmap.end(); ++itr) delete itr->second;
  cmap.clear();
  y.clear();
  p.clear();
}

/// \brief  Displays help text
///
void DisplayHelp()
{
  mexPrintf("\nmexCoherentMackenthun MATLAB wrapper the coherent phase estimator\n");
  mexPrintf("\n(C) 2009-2012 Robby McKilliam, ITR/UniSA\n");
  mexPrintf("Compiled: %s %s\n\n", __DATE__, __TIME__);
  mexPrintf("Messages are represented as %s.\n\n", STRINGIZE_VALUE_OF(MSGTYPE));
}

/// \brief  Main function for MATLAB calls.
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
  mexAtExit(Cleanup);	
	
  //construct complex gain estimator
  if (nrhs == 5){
    std::string name = mxArrayToString(prhs[0]);
    std::vector<int> D = GetIndexArray(prhs[1]); //mexPrintf("D = "); for(unsigned int i = 0; i < D.size(); i++) mexPrintf("%i, ", D[i]);
    std::vector<int> P = GetIndexArray(prhs[2]); //mexPrintf("P = "); for(unsigned int i = 0; i < P.size(); i++) mexPrintf("%i, ", P[i]);
    std::vector<complex> pilots = GetComplexArray(prhs[3]);
    int M = (int) round( *((MSGTYPE*)mxGetData(prhs[4])) );
    delete cmap[name]; //delete item of same name if it exists 
    cmap[name] = new CoherentMackenthun(D, P, pilots, M);
    return;
  }

  //construct doppler and complex gain estimator 
  else if (nrhs == 8 || nrhs == 9){
    std::string name = mxArrayToString(prhs[0]);
    std::vector<int> D = GetIndexArray(prhs[1]); //mexPrintf("D = "); for(unsigned int i = 0; i < D.size(); i++) mexPrintf("%i, ", D[i]);
    std::vector<int> P = GetIndexArray(prhs[2]); //mexPrintf("P = "); for(unsigned int i = 0; i < P.size(); i++) mexPrintf("%i, ", P[i]);
    std::vector<complex> pilots = GetComplexArray(prhs[3]);
    int M = (int) round( *((MSGTYPE*)mxGetData(prhs[4])) );
    MSGTYPE fmin = *((MSGTYPE*) mxGetData(prhs[5]));
    MSGTYPE fmax = *((MSGTYPE*) mxGetData(prhs[6]));
    MSGTYPE Ts = *((MSGTYPE*) mxGetData(prhs[7]));
    delete cmap[name]; //delete item of same name if it exists
    if (nrhs == 8){
        cmap[name] = new CoherentMackenthunWithDoppler(D, P, pilots, M, fmin, fmax, Ts);
    } else {
        MSGTYPE searchOS = *((MSGTYPE*) mxGetData(prhs[8]));
        cmap[name] = new CoherentMackenthunWithDoppler(D, P, pilots, M, fmin, fmax, Ts, searchOS);
    }
    return;
  }

  //construct doppler, doppler rate, and complex gain estimator 
  else if (nrhs == 10 || nrhs == 11){
    std::string name = mxArrayToString(prhs[0]);
    std::vector<int> D = GetIndexArray(prhs[1]); //mexPrintf("D = "); for(unsigned int i = 0; i < D.size(); i++) mexPrintf("%i, ", D[i]);
    std::vector<int> P = GetIndexArray(prhs[2]); //mexPrintf("P = "); for(unsigned int i = 0; i < P.size(); i++) mexPrintf("%i, ", P[i]);
    std::vector<complex> pilots = GetComplexArray(prhs[3]);
    int M = (int) round( *((MSGTYPE*)mxGetData(prhs[4])) );
    MSGTYPE fmin = *((MSGTYPE*) mxGetData(prhs[5]));
    MSGTYPE fmax = *((MSGTYPE*) mxGetData(prhs[6]));
    MSGTYPE frmin = *((MSGTYPE*) mxGetData(prhs[7]));
    MSGTYPE frmax = *((MSGTYPE*) mxGetData(prhs[8]));
    MSGTYPE Ts = *((MSGTYPE*) mxGetData(prhs[9]));
    delete cmap[name]; //delete item of same name if it exists
    if (nrhs == 10){
        cmap[name] = new CoherentMackenthunWithDopplerAndDopplerRate(D, P, pilots, M, fmin, fmax, frmin, frmax, Ts);
    } else {
        MSGTYPE searchOS = *((MSGTYPE*) mxGetData(prhs[10]));
        cmap[name] = new CoherentMackenthunWithDopplerAndDopplerRate(D, P, pilots, M, fmin, fmax, frmin, frmax, Ts, searchOS);
    }
    return;
  }

  //run estimator mode, pilots left to default
  else if (nrhs == 2){

    if( nlhs > 3 ) { DisplayHelp(); return; }

    //get the estimator from the map
    std::string name = mxArrayToString(prhs[0]);
    //if the estimator hasn't been created, tell the user 
    if(cmap.find(name) == cmap.end()){
      mexPrintf("Estimator named '%s' has not been created \n\n", name.c_str());
      return;
    }
    //otherwise get the estimator
    PhaseEstimator& est = *cmap[name];

    //get the signal and run the estimator
    CopyMatlabComplexToVector(prhs[1], y);
    est.estimate(y);
     
    //mexPrintf("mex out chat = %f + %f i\n", real(est.complexGainEstimate()), imag(est.complexGainEstimate()));

    //set the estimator output 
    if( nlhs > 0 ) plhs[0] = ToMatlabComplex(est.complexGainEstimate());
    if( nlhs > 1 ) plhs[1]= mxCreateDoubleScalar(est.frequencyOffsetEstimate());
    if( nlhs > 2 ) plhs[2]= mxCreateDoubleScalar(est.frequencyRateEstimate());
    
    return;
  }

  // check basic parameter format
  else {
    DisplayHelp();
    return;
  }

}
