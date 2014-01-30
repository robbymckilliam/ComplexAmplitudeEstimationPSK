/* 
 * File:   montecarlo.cpp
 * Author: Robby McKilliam
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "CoherentMackenthun.h"
#include "LDPCDec.h"
#include <random>
#include <functional>
#include <future>
#include <chrono>

using namespace std;

typedef std::chrono::duration<double, std::ratio<1>> seconds;

static constexpr unsigned int M = 2; //BPSK (you can't actually change this!)
static constexpr unsigned int numerrs = 10; //run until numerrs errors occur
static constexpr double absD = 128;

void runsim(const vector<double>& snrdbs, string name) {

  cout << "Running simulation " << name << flush;
  auto starttime = chrono::system_clock::now();

  vector<future<double>> msearray(snrdbs.size()); //store the mse

  auto endtime = chrono::system_clock::now();
  auto elapsed_seconds = std::chrono::duration_cast<seconds>(endtime-starttime).count();
  cout << " finished in  " << elapsed_seconds << " seconds." << endl;

}

int main(int argc, char** argv) {
    
    vector<double> snrdbs;//snrs in db we will run
    for(int db = -3; db <= 7; db+=1) snrdbs.push_back(db);
  
    runsim(snrdbs,"perfectchannel");
    
    
}
