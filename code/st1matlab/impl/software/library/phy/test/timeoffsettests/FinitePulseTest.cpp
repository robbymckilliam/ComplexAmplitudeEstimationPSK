/* 
 * File:   FinitePulseTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 23/01/2013, 2:19:26 PM
 */

#include "mex.h"
#include <stdlib.h>
#include <iostream>
#include "FinitePulse.h"
#include "Util.h"

using namespace TimeOffset;


bool finiteSinc4zeros() {
    std::cout << "finiteSinc with 4 zeros ... ";
    VALTYPE tol = 1e-9;
    VALTYPE nzs = 4;
    VALTYPE T = 1;
    TruncatedSincPulse sincpulse(T, nzs);
    bool pass = true;
    pass &= fabs(sincpulse.tmin() - -T * nzs) < tol;
    pass &= fabs(sincpulse.tmax() - T * nzs) < tol;
    pass &= abs(sincpulse.pulse(-4.4) - complex<VALTYPE>(0.0, 0.0)) < tol; //output from scala
    pass &= abs(sincpulse.pulse(-3.9) - complex<VALTYPE>(-0.025221324181627227, 0.0)) < tol;
    pass &= abs(sincpulse.pulse(-3.4) - complex<VALTYPE>(-0.08903843866360674, 0.0)) < tol;
    pass &= abs(sincpulse.pulse(-2.9) - complex<VALTYPE>(0.03391833252011936, 0.0)) < tol;
    pass &= abs(sincpulse.pulse(-2.4) - complex<VALTYPE>(0.12613778810677617, 0.0)) < tol;
    pass &= abs(sincpulse.pulse(-1.9) - complex<VALTYPE>(-0.05177008647807704, 0.0)) < tol;
    pass &= abs(sincpulse.pulse(-1.9) - complex<VALTYPE>(-0.05177008647807704, 0.0)) < tol;
    pass &= abs(sincpulse.pulse(0.0) - complex<VALTYPE>(1.0, 0.0)) < tol;
    pass &= abs(sincpulse.pulse(4e-3) - complex<VALTYPE>(0.9999973333354667, 0.0)) < tol;
    pass &= abs(sincpulse.pulse(1.4) - complex<VALTYPE>(-0.21623620818304484, 0.0)) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

//void normalisedTruncatedSinc() {
//    std::cout << "normalised truncated sinc ... ";
//    VALTYPE tol = 1e-4;
//    unsigned int nzs = 4;
//    VALTYPE T = 1;
//    NormalisedFinitePulse<TruncatedSincPulse> p(TruncatedSincPulse(T, nzs));
//    auto f = [&p](VALTYPE t){return std::norm(p.pulse(t));};
//    VALTYPE energy1 = trapezoidal(f, p.tmin() - 0.3, p.tmax() + 0.3, 100000); //compute pulse energy
//    VALTYPE energy2 = trapezoidal(f, p.tmin(), p.tmax(), 100000); //compute pulse energy
//    bool pass = true;
//    pass &= fabs(1.0 - energy1) < tol;
//    pass &= fabs(1.0 - energy2) < tol;
//    if (pass) std::cout << "PASS" << std::endl;
//    else std::cout << "FAIL" << std::endl;
//}

bool finiteSinc2zeros() {
    std::cout << "finiteSinc with 2 zeros ... ";
    VALTYPE tol = 1e-9;
    int nzs = 2;
    VALTYPE T = 1;
    TruncatedSincPulse sincpulse(T, nzs);
    bool pass = true;
    pass &= fabs(sincpulse.tmin() - -T * nzs) < tol;
    pass &= fabs(sincpulse.tmax() - T * nzs) < tol;
    pass &= abs(sincpulse.pulse(1.2) - complex<VALTYPE>(-0.15591488063143982, 0.0)) < tol; //output from scala
    pass &= abs(sincpulse.pulse(1.0) - complex<VALTYPE>(0.0, 0.0)) < tol; //output from scala
    pass &= abs(sincpulse.pulse(0.8) - complex<VALTYPE>(0.23387232094715982, 0.0)) < tol; //output from scala
    pass &= abs(sincpulse.pulse(0.6) - complex<VALTYPE>(0.5045511524271047, 0.0)) < tol;
    pass &= abs(sincpulse.pulse(0.4) - complex<VALTYPE>(0.756826728640657, 0.0)) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

bool durationConstructedSinc() {
    std::cout << "construct sinc pulse with duration ... ";
    VALTYPE tol = 1e-9;
    VALTYPE duration = 10.4;
    int nzs = (int) ceil(duration/2);
    VALTYPE T = 1;
    TruncatedSincPulse sincpulse(T, duration, true);
    TruncatedSincPulse sincpulsetester(T, nzs);
    bool pass = true;
    for(VALTYPE t = -duration; t < duration; t+=0.01) 
        pass &= abs(sincpulse.pulse(t) - sincpulsetester.pulse(t)) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void testRootRaisedCosineNearZero() {
    std::cout << "root raised cosine near zero... ";
    VALTYPE tol = 1e-7;
    VALTYPE nzs = 4;
    VALTYPE T = 1.0;
    VALTYPE beta = 0.5;
    TruncatedRootRaisedCosine p(T, beta, nzs);
    bool pass = true;
    pass &= fabs(1.136619772 - p.rootraisedcosine(0.0)) < tol; //from mathematica
    pass &= fabs(1.136617045 - p.rootraisedcosine(1e-3)) < tol; //from mathematica
    pass &= fabs(1.136617045 - p.rootraisedcosine(-1e-3)) < tol; //from mathematica
    pass &= fabs(1.136576129 - p.rootraisedcosine(-4e-3)) < tol; //from mathematica
    pass &= fabs(1.136576129 - p.rootraisedcosine(4e-3)) < tol; //from mathematica
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void testRootRaisedCosineNearBeta4() {
    std::cout << "root raised cosine near beta/4... ";
    VALTYPE tol = 1e-6;
    VALTYPE nzs = 4;
    VALTYPE T = 1.0;
    VALTYPE beta = 0.5;
    VALTYPE t = 1.0 / 4 / beta;
    TruncatedRootRaisedCosine p(T, beta, nzs);
    bool pass = true;
    pass &= fabs(0.5786324696 - p.rootraisedcosine(t)) < tol; //from Mathematica
    pass &= fabs(0.5786324696 - p.rootraisedcosine(-t)) < tol; //from Mathematica

    pass &= fabs(0.5784538720 - p.rootraisedcosine(t + 1e-4)) < tol; //from Mathematica
    pass &= fabs(0.5788110635 - p.rootraisedcosine(t - 1e-4)) < tol; //from Mathematica
    pass &= fabs(0.5784538720 - p.rootraisedcosine(-t - 1e-4)) < tol; //from Mathematica
    pass &= fabs(0.5788110635 - p.rootraisedcosine(-t + 1e-4)) < tol; //from Mathematica

    pass &= fabs(0.5793468225 - p.rootraisedcosine(t - 4e-4)) < tol; //from Mathematica
    pass &= fabs(0.5779180564 - p.rootraisedcosine(t + 4e-4)) < tol; //from Mathematica
    pass &= fabs(0.5793468225 - p.rootraisedcosine(-t + 4e-4)) < tol; //from Mathematica
    pass &= fabs(0.5779180564 - p.rootraisedcosine(-t - 4e-4)) < tol; //from Mathematica
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void testRootRaisedCosine() {
    std::cout << "root raised cosine... ";
    VALTYPE tol = 1e-6;
    VALTYPE nzs = 4;
    VALTYPE T = 1.0;
    VALTYPE beta = 0.5;
    TruncatedRootRaisedCosine p(T, beta, nzs);
    VALTYPE d[9] = {2.0 / 15 / pi, -1.0 / 3 / sqrt(2.0) / pi, -1.0 / 3 / pi, (pi + 2) / 2 / sqrt(2.0) / pi, 0.5 + 2 / pi, (pi + 2) / 2 / sqrt(2.0) / pi, -1.0 / 3 / pi, -1.0 / 3 / sqrt(2.0) / pi, 2.0 / 15 / pi}; //from Mathematica
    std::vector<VALTYPE> expected ;
    for(int i = 0; i < 9; i++) expected.push_back(d[i]);
    std::vector<VALTYPE> test;
    for (VALTYPE t = -2.0; t <= 2.0; t += 0.5) test.push_back(p.rootraisedcosine(t));
    bool pass = true;
    for (int i = 0; i < test.size(); i++)
        pass &= fabs(test[i] - expected[i]) < tol; //from Mathematica
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void testRootRaisedCosineFindZeros() {
    std::cout << "root raised cosine find zeros ... ";
    const VALTYPE tol = 1e-4;
    const VALTYPE T = 1.0;
    const VALTYPE beta = 0.5;
    //test against roots found by Mathematica
    bool pass = true;
    pass &= fabs(0.873584 - TruncatedRootRaisedCosine(T, beta, 1).tmax()) < tol;
    pass &= fabs(1.69555 - TruncatedRootRaisedCosine(T, beta, 2).tmax()) < tol;
    pass &= fabs(-1.69555 - TruncatedRootRaisedCosine(T, beta, 2).tmin()) < tol;
    pass &= fabs(2.35734 - TruncatedRootRaisedCosine(T, beta, 3).tmax()) < tol;
    pass &= fabs(2.96409 - TruncatedRootRaisedCosine(T, beta, 4).tmax()) < tol;
    pass &= fabs(3.68054 - TruncatedRootRaisedCosine(T, beta, 5).tmax()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

bool durationConstructedRRC() {
    std::cout << "construct rrc pulse with duration ... ";
    VALTYPE tol = 1e-9;
    const VALTYPE beta = 1.0/3.0;
    VALTYPE duration = 10;
    int nzs = 7; //obtained by inspection.
    VALTYPE T = 1;
    TruncatedRootRaisedCosine rrcpulse(T, beta, duration, true);
    TruncatedRootRaisedCosine rrcpulsetester(T, beta, nzs);
    bool pass = true;
    for(VALTYPE t = -duration; t < duration; t+=0.01) 
        pass &= abs(rrcpulse.pulse(t) - rrcpulsetester.pulse(t)) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void testWithLookUpTable() {
    std::cout << "look up table pulse ... ";
    const VALTYPE tol = 1e-3;
    const int nzs = 4;
    const VALTYPE T = 1.0;
    const VALTYPE beta = 0.5;
    TruncatedRootRaisedCosine tx(T,beta,nzs);
    FinitePulseWithLookupTable<TruncatedRootRaisedCosine> luttx(tx,100000);
    bool pass = true;
    for( VALTYPE t = tx.tmin() - 1.0; t < tx.tmax() + 1.0; t+=0.01)
      pass &= std::abs<VALTYPE>(tx.pulse(t) - luttx.pulse(t)) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
  }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    finiteSinc4zeros();
    finiteSinc2zeros();
    //normalisedTruncatedSinc();
    durationConstructedSinc();

    testRootRaisedCosineNearZero();
    testRootRaisedCosineNearBeta4();
    testRootRaisedCosine();
    testRootRaisedCosineFindZeros();
    testWithLookUpTable();
    durationConstructedRRC();
}

