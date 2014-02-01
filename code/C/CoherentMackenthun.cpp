#include "CoherentMackenthun.h"
#include <math.h>
#include <float.h>
#include <algorithm>
#include <sstream> 
#include <iostream>

#ifdef _WIN32
#define log2(x) 			(log(x)/log((double)2.0))
#define round(x)			(floor((x)+0.5))
#endif

/// \brief	Standard constructor.
CoherentMackenthun::CoherentMackenthun(
 				       const std::vector<int>& Din, 
 				       const std::vector<int>& Pin,
 				       const std::vector<complexd>& pin,
  				       const unsigned int Min) :
  L(Din.size() + Pin.size()),
  M(Min),
  absD(Din.size()),
  absP(Pin.size()),
  D(Din),
  P(Pin),
  p(pin),
  w(2*pi/M),
  eta(std::polar<double>(1.0,w)),
  nu(std::polar<double>(1.0,w) - complexd(1.0,0.0))
{
  if( pin.size() != P.size() ) throw std::string("The number of pilot symbols does not match the number of indices in Pin");
  if( M < 1 )  throw std::string("M must be a positive integer");
  //allocate memory
  z.resize(D.size(), IndexedReal(0,0));
  for(unsigned int i = 0; i < z.size(); i++) z[i] = IndexedReal(0,i); //initialise z with the identity permutation 
  g.resize(D.size(), complexd(0,0));
  PUD.resize(L, 0); //union P and D
  for(unsigned int i = 0; i < P.size(); i++) PUD[i] = P[i];
  for(unsigned int i = 0; i < D.size(); i++) PUD[i + P.size()] = D[i];
  max_index = *std::max_element(PUD.begin(), PUD.end());
  min_index = *std::min_element(PUD.begin(), PUD.end());
}

/// Run the estimator on data y
void CoherentMackenthun::estimate(const std::vector<complexd>& y)
{
  estimate(y, p);
}

/// Run the estimator with data y and pilot symbols p
void CoherentMackenthun::estimate(const std::vector<complexd>& y, const std::vector<complexd>& p)
{
  if( y.size() <= max_index ) {
    std::ostringstream str;
    str << "Exception in function CoherentMackenthun::estimate. The data length y is ";
    str << y.size() << " but the maximum symbol index is " << max_index;
    throw str.str();
  }
  if( p.size() != P.size() ) {
    std::ostringstream str;
    str << "Exception in function CoherentMackenthun::estimate. The number of pilot symbols given is ";
    str << p.size() << " but the number of indices is " << P.size();
    throw str.str();
  }
  
  //the value A, norm of received signal
  double A = 0.0;
  for(int i = 0; i < L; i++) A += std::norm(y[PUD[i]]);
  
  //setup sequences and sort
  complexd Y(0.0,0.0);
  for(unsigned int i=0; i < P.size(); i++) Y += y[P[i]]*conj(p[i]);
  for(unsigned int i=0; i < D.size(); i++){
    double phi = std::arg(y[D[i]]);
    double u = w*round(phi/w);
    g[i] = y[D[i]] * std::polar<double>(1,-u);
    z[i] = IndexedReal(phi-u, i);
    Y += g[i];
  }
  std::sort(z.begin(), z.end());

  chat = Y / ((double) L);
  Qhat = std::norm(Y) / L;
  noisevarhat = (A - Qhat) / L;
  for(unsigned int k = 0; k < M*D.size(); k++ ){
    Y += nu*g[sigma(k)];
    g[sigma(k)] *= eta;
    double Q = std::norm(Y) / L;
    if( Q > Qhat ){
      chat = Y / ((double) L);
      Qhat = Q;
      noisevarhat = (A - Qhat) / L;
    }
  }
  
}