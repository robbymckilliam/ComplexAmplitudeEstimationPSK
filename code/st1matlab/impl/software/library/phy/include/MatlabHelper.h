//
//  MatlabHelper.h
//  QAMModulator
//
//  Created by Gottfried Lechner on 12/04/12.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#ifndef _CoherehentMackenthun_MatlabHelper_h
#define _CoherehentMackenthun_MatlabHelper_h

#ifndef MSGTYPE
  #define MSGTYPE double
#endif

#include <complex>
#include <vector>
#include "mex.h"

typedef std::complex<MSGTYPE> complexT;

std::vector<MSGTYPE> GetFloatArray(const mxArray *src);
std::vector<complexT> GetComplexArray(const mxArray *src);
std::vector<int> GetIntArray(const mxArray *src);
std::vector<int> GetIndexArray(const mxArray *src);
mxArray* ToMatlabComplex(complexT src);
mxArray* ToMatlabComplex(std::vector<complexT> src);
void CopyMatlabComplexToVector(const mxArray *src, std::vector<complexT>& v);
void CopyMatlabComplexToVector(const mxArray *src, std::vector<MSGTYPE>& v);
void CopyMatlabComplexToVector(const mxArray *src, std::vector<int>& v);

void GetFloatArray(const mxArray *src, MSGTYPE *dst, bool imag = false);
void GetIntArray(const mxArray *src, unsigned int *dst);
void SetFloatArray(MSGTYPE *src, mxArray *dst, bool imag = false);
void SetIntArray(unsigned int *src, mxArray *dst);

#endif
