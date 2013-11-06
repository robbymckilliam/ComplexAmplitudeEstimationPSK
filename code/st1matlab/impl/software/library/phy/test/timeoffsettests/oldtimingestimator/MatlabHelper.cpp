//
//  MatlabHelper.cpp
//  CoherentMackenthun
//
//  Created by Gottfried Lechner on 12/04/12.
//  Modified by Robby McKilliam on 5/04/12
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#include "MatlabHelper.h"

#ifdef _WIN32
#define round(x)			(floor((x)+0.5))
#endif

/// \brief  Copy a matlab float array into a vector of MSGTYPE (i.e. doubles or floats).  Ignores imaginary part.
///
std::vector<MSGTYPE> GetFloatArray(const mxArray *src)
{
  unsigned int N = mxGetNumberOfElements(src);
  std::vector<MSGTYPE> ret(N);
  
  double  *dr;
  float   *fr;
  
  switch(mxGetClassID(src))
  {
    case mxDOUBLE_CLASS:
      dr = (double*) mxGetData(src);
      for(unsigned int n = 0; n < N; n++) ret[n] = (MSGTYPE)dr[n];
      break;
    case mxSINGLE_CLASS:
      fr = (float*) mxGetData(src);
      for(unsigned int n = 0; n < N; n++) ret[n] = (MSGTYPE)fr[n];
      break;
  }
  return ret;
}

//Copies the values in a matlab mxArray* into a std::vector<complex>.  Will resize the std::vector
//if neccessary. 
void CopyMatlabComplexToVector(const mxArray *src, std::vector<complex>& v){
  unsigned int N = mxGetNumberOfElements(src);
  if( v.size() != N ) v.resize(N, complex(0,0));
  double* r = (double*) mxGetData(src);
  double* i = (double*) mxGetImagData(src);
  if(i != NULL) for(unsigned int n = 0; n < N; n++) v[n] = complex(r[n],i[n]); 
  else for(unsigned int n = 0; n < N; n++) v[n] = complex(r[n],0.0);
}

//Copies the real part of a matlab mxArray* into a std::vector<MSGTYPE>.  Will resize the std::vector
//if neccessary. 
void CopyMatlabComplexToVector(const mxArray *src, std::vector<MSGTYPE>& v){
  unsigned int N = mxGetNumberOfElements(src);
  if( v.size() != N ) v.resize(N, 0.0);
  double* r = (double*) mxGetData(src); 
  for(unsigned int n = 0; n < N; n++) v[n] = r[n];
}

//Copies the real part of a matlab mxArray* into a std::vector<int>.  Will resize the std::vector
//if neccessary. 
void CopyMatlabComplexToVector(const mxArray *src, std::vector<int>& v){
  unsigned int N = mxGetNumberOfElements(src);
  if( v.size() != N ) v.resize(N, 0);
  double* r = (double*) mxGetData(src); 
  for(unsigned int n = 0; n < N; n++) v[n] = (int) round(r[n]);
}

/// \brief  Copy a complex matlab array into a vector of std::complex.
///
std::vector<complex> GetComplexArray(const mxArray *src)
{
  unsigned int N = mxGetNumberOfElements(src);
  std::vector<complex> ret(N);
  
  double  *dr, *di;
  float   *fr, *fi;
  
  switch(mxGetClassID(src))
  {
    case mxDOUBLE_CLASS:
      dr = (double*) mxGetData(src);
      di = (double*) mxGetImagData(src);
      if(di != NULL) for(unsigned int n = 0; n < N; n++) ret[n] = complex(dr[n],di[n]); 
      else  for(unsigned int n = 0; n < N; n++) ret[n] = complex(dr[n],0.0);
      break;
    case mxSINGLE_CLASS:
      fr = (float*) mxGetData(src);
      fi = (float*) mxGetImagData(src);
      if(fi != NULL) for(unsigned int n = 0; n < N; n++) ret[n] = complex(fr[n],fi[n]); 
      else  for(unsigned int n = 0; n < N; n++) ret[n] = complex(fr[n],0.0);
      break;
  }
      
  return ret;
}

/// \brief  Copy a complex to a matlab complex.
///
mxArray* ToMatlabComplex(std::complex<MSGTYPE> src)
{
  mxArray *p = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
  MSGTYPE *ar = (MSGTYPE*) mxGetData(p); 
  ar[0] = real(src);
  MSGTYPE *ai = (MSGTYPE*) mxGetImagData(p); 
  ai[0] = imag(src);
  return p;
}

/// \brief  Copy a complex vector to a matlab complex vector.
///
mxArray* ToMatlabComplex(std::vector<std::complex<MSGTYPE> > src)
{
  mxArray *p = mxCreateDoubleMatrix(1, src.size(), mxCOMPLEX);
  double  *ar = (double*) mxGetData(p); 
  double  *ai = (double*) mxGetImagData(p);
  for(unsigned int n=0;n<src.size();n++)
  {
    ar[n] = src[n].real();
    ai[n] = src[n].imag();
  }
  return p;
}


/// \brief  Copy a matlab float (or complex) array into a vector of ints.  Ignores imaginary part.
///
std::vector<int> GetIntArray(const mxArray *src)
{
  int N = mxGetNumberOfElements(src);
  std::vector<int > ret(N);
  MSGTYPE* r = (MSGTYPE*) mxGetData(src);
  for(int n = 0; n < N; n++) ret[n] = (int)round(r[n]); 
  return ret;
}

/// \brief  Copy a matlab float (or complex) array into a vector of ints but subtracts one from each.  
/// This maps the indices in matlab 1,2,3,.... to C indices 0,1,2,... . Ignores imaginary part.
///
std::vector<int> GetIndexArray(const mxArray *src)
{
  int N = mxGetNumberOfElements(src);
  std::vector<int > ret(N);
  MSGTYPE* r = (MSGTYPE*) mxGetData(src);
  for(int n = 0; n < N; n++) ret[n] = ((int)round(r[n])) - 1; 
  return ret;
}

/// \brief  Performs copy&convert operation on floating point input data.
///
void GetFloatArray(const mxArray *src, MSGTYPE *dst, bool imag)
{
  double        *dptr;
  float         *fptr;
  unsigned int  n;
  
  switch(mxGetClassID(src))
  {
    case mxDOUBLE_CLASS:
      if(imag)
        dptr = (double*)mxGetImagData(src);
      else
        dptr = (double*)mxGetData(src);
      for(n=0;n<mxGetNumberOfElements(src);n++)
        dst[n] = (MSGTYPE)dptr[n];
      break;
    case mxSINGLE_CLASS:
      if(imag)
        fptr = (float*)mxGetImagData(src);
      else
        fptr = (float*)mxGetData(src);
      for(n=0;n<mxGetNumberOfElements(src);n++)
        dst[n] = (MSGTYPE)fptr[n];
      break;
    default:
      mexErrMsgTxt("Wrong input data type. Only single and double are supported.");
  }
}


/// \brief  Performs copy&convert operation on integer input data.
///
void GetIntArray(const mxArray *src, unsigned int *dst)
{
  double        *dptr;
  float         *fptr;
  unsigned int  n;
  
  switch(mxGetClassID(src))
  {
    case mxDOUBLE_CLASS:
      dptr = (double*)mxGetData(src);
      for(n=0;n<mxGetNumberOfElements(src);n++)
        dst[n] = (int)dptr[n];
      break;
    case mxSINGLE_CLASS:
      fptr = (float*)mxGetData(src);
      for(n=0;n<mxGetNumberOfElements(src);n++)
        dst[n] = (int)fptr[n];
      break;
    default:
      mexErrMsgTxt("Wrong input data type. Only single and double are supported.");
  }
}

/// \brief  Performs copy&convert operation on floating point output data.
///
void SetFloatArray(MSGTYPE *src, mxArray *dst, bool imag)
{
  double        *dptr;
  float         *fptr;
  unsigned int  n;
  
  switch(mxGetClassID(dst))
  {
    case mxDOUBLE_CLASS:
      if(imag)
        dptr = (double*)mxGetImagData(dst);
      else
        dptr = (double*)mxGetData(dst);
      for(n=0;n<mxGetNumberOfElements(dst);n++)
        dptr[n] = (double)src[n];
      break;
    case mxSINGLE_CLASS:
      if(imag)
        fptr = (float*)mxGetImagData(dst);
      else
        fptr = (float*)mxGetData(dst);
      for(n=0;n<mxGetNumberOfElements(dst);n++)
        fptr[n] = (float)src[n];
      break;
    default:
      mexErrMsgTxt("Wrong input data type. Only single and double are supported.");
  }
}

/// \brief  Performs copy&convert operation on integer output data.
///
void SetIntArray(unsigned int *src, mxArray *dst)
{
  double        *dptr;
  float         *fptr;
  unsigned int  n;
  
  switch(mxGetClassID(dst))
  {
    case mxDOUBLE_CLASS:
      dptr = (double*)mxGetData(dst);
      for(n=0;n<mxGetNumberOfElements(dst);n++)
        dptr[n] = (double)src[n];
      break;
    case mxSINGLE_CLASS:
      fptr = (float*)mxGetData(dst);
      for(n=0;n<mxGetNumberOfElements(dst);n++)
        fptr[n] = (float)src[n];
      break;
    default:
      mexErrMsgTxt("Wrong input data type. Only single and double are supported.");
  }
}
