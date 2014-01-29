#include "mex.h"
#include "CoherentMackenthun.h"
#include "MatlabHelper.h"
#include <math.h>

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

static const double TOL = 0.0000001;

/// \brief  Main function for MATLAB calls.
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
  
  bool pass;

  mexPrintf("Testing GetFloatArray ... ");
  mwSize N = 5;
  mxArray *p = mxCreateDoubleMatrix(N,1,mxREAL);
  MSGTYPE *d = (MSGTYPE*) mxGetData(p);
  for(int n = 0; n < N; n++) d[n] = (MSGTYPE) n;
  std::vector<MSGTYPE> v1 = GetFloatArray(p); 
  pass = true;
  for(int n = 0; n < N; n++) pass &= abs(v1[n] - d[n]) < TOL;
  if(pass) mexPrintf("PASS\n");
  else mexPrintf("FAIL\n");

  mexPrintf("Testing GetComplexArray ... ");
  p = mxCreateDoubleMatrix(N,1,mxCOMPLEX);
  MSGTYPE *dr = (MSGTYPE*) mxGetData(p);
  MSGTYPE *di = (MSGTYPE*) mxGetImagData(p);
  for(int n = 0; n < N; n++) {
    dr[n] = (MSGTYPE) n;
    di[n] = (MSGTYPE) -10.1*n;
  }
  std::vector<std::complex<MSGTYPE> > v2 = GetComplexArray(p); 
  pass = true;
  for(int n = 0; n < N; n++) pass &= abs(v2[n] - std::complex<MSGTYPE>(dr[n], di[n])) < TOL;
  if(pass) mexPrintf("PASS\n");
  else mexPrintf("FAIL\n");
				
  mexPrintf("Testing GetIntArray ... ");
  p = mxCreateDoubleMatrix(N,1,mxREAL);
  d = (MSGTYPE*) mxGetData(p);
  for(int n = 0; n < N; n++) d[n] = (MSGTYPE) n;
  std::vector<int> v3 = GetIntArray(p); 
  pass = true;
  for(int n = 0; n < N; n++) pass &= abs(v3[n] - (int)d[n]) < TOL;
  if(pass) mexPrintf("PASS\n");
  else mexPrintf("FAIL\n");	 

  mexPrintf("Testing GetIndexArray ... ");
  p = mxCreateDoubleMatrix(N,1,mxREAL);
  d = (MSGTYPE*) mxGetData(p);
  for(int n = 0; n < N; n++) d[n] = (MSGTYPE) n;
  std::vector<int> v4 = GetIndexArray(p); 
  pass = true;
  for(int n = 0; n < N; n++) pass &= (v4[n] - (int)d[n] + 1) == 0;
  //for(int n = 0; n < N; n++) mexPrintf("vals %i %i \n", v4[n], (int)d[n]);
  if(pass) mexPrintf("PASS\n");
  else mexPrintf("FAIL\n");	

  mexPrintf("Testing CopyMatlabComplexToVector (complex) ... ");
  p = mxCreateDoubleMatrix(N,1,mxCOMPLEX);
  dr = (MSGTYPE*) mxGetData(p);
  di = (MSGTYPE*) mxGetImagData(p);
  for(int n = 0; n < N; n++) {
    dr[n] = (MSGTYPE) n;
    di[n] = (MSGTYPE) -10.1*n;
  }
  std::vector<complex> v5;
  CopyMatlabComplexToVector(p,v5);
  pass = true;
  for(int n = 0; n < N; n++) pass &= abs(v5[n] - std::complex<MSGTYPE>(dr[n], di[n])) < TOL;
  //for(int n = 0; n < N; n++) mexPrintf("(%f, %f), (%f, %f)", std::real<MSGTYPE>(v5[n]), std::imag<MSGTYPE>(v5[n]), dr[n], di[n]);
  if(pass) mexPrintf("PASS\n");
  else mexPrintf("FAIL\n");

mexPrintf("Testing CopyMatlabComplexToVector (float) ... ");
  p = mxCreateDoubleMatrix(N,1,mxCOMPLEX);
  dr = (MSGTYPE*) mxGetData(p);
  di = (MSGTYPE*) mxGetImagData(p);
  for(int n = 0; n < N; n++) {
    dr[n] = (MSGTYPE) n;
    di[n] = (MSGTYPE) -10.1*n;
  }
  std::vector<double> v6;
  CopyMatlabComplexToVector(p,v6);
  pass = true;
  for(int n = 0; n < N; n++) pass &= abs(v6[n] - dr[n]) < TOL;
  if(pass) mexPrintf("PASS\n");
  else mexPrintf("FAIL\n");

mexPrintf("Testing CopyMatlabComplexToVector (int) ... ");
  p = mxCreateDoubleMatrix(N,1,mxCOMPLEX);
  dr = (MSGTYPE*) mxGetData(p);
  di = (MSGTYPE*) mxGetImagData(p);
  for(int n = 0; n < N; n++) {
    dr[n] = (MSGTYPE) n;
    di[n] = (MSGTYPE) -10.1*n;
  }
  std::vector<int> v7;
  CopyMatlabComplexToVector(p,v7);
  pass = true;
  for(int n = 0; n < N; n++) pass &= abs(v7[n] - dr[n]) < TOL;
  if(pass) mexPrintf("PASS\n");
  else mexPrintf("FAIL\n");

  mexPrintf("Testing toMatlabComplex ... ");
  std::complex<MSGTYPE> c(1.0,-2.3);
  p = ToMatlabComplex(c);
  dr = (MSGTYPE*) mxGetData(p);
  di = (MSGTYPE*) mxGetImagData(p);
  pass &= abs( dr[0] - real(c) ) < TOL;
  pass &= abs( di[0] - imag(c) ) < TOL;
  //mexPrintf("vals %f %f %f %f \n", dr[0], di[0], real(c), imag(c));
  if(pass) mexPrintf("PASS\n");
  else mexPrintf("FAIL\n");

  mexPrintf("Testing Coherent Mackenthun ... ");
  int M = 4; //QPSK
  std::vector<int> P;
  P.push_back(0);
  std::vector<int> D;
  D.push_back(1); D.push_back(2); D.push_back(3);  D.push_back(4);
  std::vector<complex> pilots;
  pilots.push_back(complex(1,0));
  std::vector<complex> y; 
  y.push_back(complex(-0.15662206555348682, 0.998001546268412));
  y.push_back(complex(-0.15571514521240135, 0.993262721949701)); 
  y.push_back(complex(0.9856798266529281, 0.12784677325244354));
  y.push_back(complex(0.14524463440340715, -0.9886039113420713));
  y.push_back(complex(0.14550394304366623, -0.9934311310398819));
  CoherentMackenthun cmack(D, P, pilots, M); 
  cmack.estimate(y);
  complex ahat = cmack.complexGainEstimate();
  complex e1(-0.14618651229308086, 0.9917958274505988);
  if(abs(e1 - ahat) < TOL) mexPrintf("PASS\n");
  else mexPrintf("FAIL\n");

  mexPrintf("\n");

}
