/* 
 * File:   BrentTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 26/01/2013, 3:48:58 PM
 */

#include "mex.h"
#include <stdlib.h>
#include <iostream>
#include "Util.h"

using namespace TimeOffset;

VALTYPE quadfunc(VALTYPE x) { return (x-2.0)*(x-2.0); } 
  void QuadraticTest() {
    VALTYPE tol = 1e-7;
   Brent<VALTYPE(VALTYPE)> opt(quadfunc, -4.0, 1.0, 5.0);
    std::cout << "Brent QuadraticTest ... ";
    bool pass = fabs(0.0 - opt.fmin()) < tol;
    pass &= fabs(2.0 - opt.xmin()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL " << opt.fmin() << " " << opt.xmin() << std::endl;
  }
  
  VALTYPE quartfunc(VALTYPE x)  { return x*x*x*x; }
    void QuarticTest() {
    VALTYPE tol = 1e-7;
    Brent<VALTYPE(VALTYPE)> opt(quartfunc, -2.0, 1.0, 2.0);
    std::cout << "Brent QuarticTest ... ";
    bool pass = fabs(0.0 - opt.fmin()) < tol;
    pass &= fabs(0.0 - opt.xmin()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL " << opt.fmin() << " " << opt.xmin() << std::endl;
  }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    QuadraticTest();
    QuarticTest();
}

