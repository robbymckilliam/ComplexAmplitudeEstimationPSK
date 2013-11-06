/* 
 * File:   UtilTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 24/01/2013, 11:21:20 AM
 */

#include "mex.h"
#include <stdlib.h>
#include <iostream>
#include "Util.h"

using namespace TimeOffset;


VALTYPE sqrint(VALTYPE a, VALTYPE b) { return b*b*b/3 - a*a*a/3; }
void trapezoidalSquareTest() {
    VALTYPE tol = 1e-6;
    std::cout << "Trapezoidal Square Test ... ";
    unsigned int N = 5000;
    bool pass = abs(trapezoidal(sqr, -3.0, 4.0, N) - sqrint(-3.0,4.0)) < tol;   
    if(pass) std::cout << "PASS" << std::endl;
    else {
        std::cout << "FAIL" << std::endl;
        std::cout << "expected " << sqrint(-3.0,4.0) << ", but was " << trapezoidal(sqr, -3.0, 4.0, N) << std::endl;
    }
}

void gcdTest() {
    std::cout << "gcd ... ";
    bool pass = gcd(10,2)==2;
    pass &= gcd(3458,4864)==38;
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void extended_gcd_Test() {
    std::cout << "extended gcd ... ";
    int gcd, x, y;
    bool pass = true;
    extended_gcd(10,2,gcd,x,y);
    pass &= (gcd==2)&&((10*x + 2*y) == gcd);
    extended_gcd(3458,4864,gcd,x,y);
    pass &= (gcd==38)&&((3458*x + 4864*y) == gcd);
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

VALTYPE feqx(VALTYPE x) {return x;}
void fzeroLinearTest() {
    std::cout << "Bisection linear ... ";
    VALTYPE tol = 1e-7;
    VALTYPE x = Bisection<VALTYPE(VALTYPE)>(feqx, -11.0,8.0,tol).zero();
    bool pass = fabs(0.0 - feqx(x)) <  tol;
    pass &= fabs(0.0 - x) <  tol;
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
  }

VALTYPE cubfunc(VALTYPE x) {return x*x*(x-1); } 
void fzeroCubicTest() {
    std::cout << "Bisection cubic ... ";
    VALTYPE tol = 1e-7;
    VALTYPE x = Bisection<VALTYPE(VALTYPE)>(cubfunc, 0.5,1.7,tol).zero();
    bool pass = true;
    pass &= fabs(1.0 - x) <  tol;
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
  }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    trapezoidalSquareTest();
    extended_gcd_Test();
    gcdTest();
    fzeroLinearTest();
    fzeroCubicTest();
    
}

