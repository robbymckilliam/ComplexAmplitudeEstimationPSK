/* 
 * File:   newsimpletest.cpp
 * Author: Robby McKilliam
 *
 * Created on 29/01/2014, 1:38:53 PM
 */

#include <stdlib.h>
#include <iostream>
#include <functional>
#include "CoherentMackenthun.h"

#define	TOL 1e-7

using namespace std;

/*
 * Test Coherent Mackenthun code.
 */

bool testEstimate() {
    int M = 4; //QPSK
    std::vector<int> P;
    P.push_back(0);
    std::vector<int> D;
    D.push_back(1);
    D.push_back(2);
    D.push_back(3);
    D.push_back(4);
    std::vector<complexd> pilots;
    pilots.push_back(complexd(1, 0));
    std::vector<complexd> y;
    y.push_back(complexd(-0.15662206555348682, 0.998001546268412));
    y.push_back(complexd(-0.15571514521240135, 0.993262721949701));
    y.push_back(complexd(0.9856798266529281, 0.12784677325244354));
    y.push_back(complexd(0.14524463440340715, -0.9886039113420713));
    y.push_back(complexd(0.14550394304366623, -0.9934311310398819));
    CoherentMackenthun cmack(D, P, pilots, M);
    cmack.estimate(y);
    complexd ahat = cmack.complexGainEstimate();
    complexd e1(-0.14618651229308086, 0.9917958274505988);
    if (abs(e1 - ahat) < TOL) return true;
    else return false;
}

bool testEstimateNoiseVar() {
    int M = 2;
    int absP = 100;
    int absD = 50;
    int L = absP + absD;
    std::vector<int> P;
    for(int i = 0; i < absP; i++) P.push_back(i); //pilots at the front
    std::vector<int> D;
    for(int i = absP; i < absD; i++) D.push_back(i); //data at the back
}

void runtest(string name, function<bool()> test) {
    cout << name << " ... ";
    if (!test()) cout << "FAIL" << endl;
    else cout << "pass" << endl;
}

int main(int argc, char** argv) {
    runtest("test coherent Mackenthun estimate", testEstimate);
    return (EXIT_SUCCESS);
}

