#include "CoherentMackenthun.h"
//#include "mex.h"
#include <math.h>
#include <float.h>
#include <algorithm>
#include <sstream> 
#include <iostream>

#define PI 3.1415926535897932384626433

#ifdef _WIN32
#define log2(x) 			(log(x)/log((MSGTYPE)2.0))
#define round(x)			(floor((x)+0.5))
#endif

/// \brief	Standard constructor.
CoherentMackenthun::CoherentMackenthun(
 				       const std::vector<int>& Din, 
 				       const std::vector<int>& Pin,
 				       const std::vector<complex>& pin,
  				       const unsigned int Min) :
  L(Din.size() + Pin.size()),
  M(Min),
  absD(Din.size()),
  absP(Pin.size()),
  D(Din),
  P(Pin),
  p(pin),
  w(2*PI/M),
  eta(std::polar<MSGTYPE>(1.0,w)),
  nu(std::polar<MSGTYPE>(1.0,w) - complex(1.0,0.0))
{
  if( pin.size() != P.size() ) throw std::string("The number of pilot symbols does not match the number of indices in Pin");
  if( M < 1 )  throw std::string("M must be a positive integer");
  //allocate memory
  z.resize(D.size(), IndexedReal(0,0));
  for(unsigned int i = 0; i < z.size(); i++) z[i] = IndexedReal(0,i); //initialise z with the identity permutation 
  g.resize(D.size(), complex(0,0));
  PUD.resize(L, 0); //union P and D
  for(unsigned int i = 0; i < P.size(); i++) PUD[i] = P[i];
  for(unsigned int i = 0; i < D.size(); i++) PUD[i + P.size()] = D[i];
  max_index = *std::max_element(PUD.begin(), PUD.end());
  min_index = *std::min_element(PUD.begin(), PUD.end());
}

/// Run the estimator on data y
void CoherentMackenthun::estimate(const std::vector<complex>& y)
{
  estimate(y, p);
}

/// Run the estimator with data y and pilot symbols p
void CoherentMackenthun::estimate(const std::vector<complex>& y, const std::vector<complex>& p)
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

  //setup sequences and sort
  complex Y(0.0,0.0);
  for(unsigned int i=0; i < P.size(); i++) Y += y[P[i]]*conj(p[i]);
  for(unsigned int i=0; i < D.size(); i++){
    MSGTYPE phi = std::arg(y[D[i]]);
    MSGTYPE u = w*round(phi/w);
    g[i] = y[D[i]] * std::polar<MSGTYPE>(1,-u);
    z[i] = IndexedReal(phi-u, i);
    Y += g[i];
  }
  sort(z.begin(), z.end());

  MSGTYPE fL = (MSGTYPE) L;
  chat = Y / fL;
  Qhat = norm(Y) / fL;
  for(unsigned int k = 0; k < M*D.size(); k++ ){
    Y += nu*g[sigma(k)];
    g[sigma(k)] *= eta;
    MSGTYPE Q = norm(Y) / fL;
    if( Q > Qhat ){
      chat = Y / fL;
      Qhat = Q;
    }
  }
  
}

CoherentMackenthunWithDoppler::CoherentMackenthunWithDoppler(
							     const std::vector<int>& D,
							     const std::vector<int>& P,
							     const std::vector<complex>& p,
							     const unsigned int M,
							     const MSGTYPE fmin_,
							     const MSGTYPE fmax_,
							     const MSGTYPE T,
							     const MSGTYPE search_oversample) 
  : CoherentMackenthun(D,P,p,M), 
    Ts(T),
    fmin(fmin_),
    fmax(fmax_)
{
  int indif = max_index-min_index;
  int ftrials = (int)ceil(search_oversample*M*(T*indif)*(fmax - fmin));
  if( ftrials == 0 ) ftrials = 1; //atleast one trial
  fstep = (fmax - fmin)/ftrials;
  if( fstep == 0.0) fstep = 1.0; //prevent potential infinite loop
  //allocate working memory
  yf.resize(max_index+1, complex(0,0));
  X.resize(L, complex(0,0));
}

/// Run the estimator
void CoherentMackenthunWithDoppler::estimate(const std::vector<complex>& y, const std::vector<complex>& p)
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
  
  MSGTYPE Qhat = -1.0;
  chat = complex(0,0);
  for(MSGTYPE f = fmin; f <= fmax; f += fstep){

    //setup sequences and sort
    complex Y(0.0,0.0);
    for(unsigned int i=0; i < P.size(); i++) Y += y[P[i]] * conj(p[i]) * std::polar<MSGTYPE>(1.0, -2*PI*f*Ts*(P[i]+1));
    for(unsigned int i=0; i < D.size(); i++){
      int n = z[i].i; //use last permutation order (should make sort faster)
      //int n = i;
      complex DCtemp = std::polar<MSGTYPE>(1.0, -2*PI*f*Ts*(D[n]+1)); //Use D[i]+1 here and above to be consistent with MATLAB indexing
      MSGTYPE phi = std::arg( y[D[n]] * DCtemp ); 
      MSGTYPE u = w*round(phi/w);
      g[n] = y[D[n]] * std::polar<MSGTYPE>(1,-u) * DCtemp;
      z[i] = IndexedReal(phi-u, n);
      Y += g[n];
    }
    sort(z.begin(), z.end());

    MSGTYPE fL = (MSGTYPE) L;
    for(unsigned int k = 0; k < M*D.size(); k++ ){
      Y += nu*g[sigma(k)];
      g[sigma(k)] *= eta;
      MSGTYPE Q = norm(Y) / fL;
      if( Q > Qhat ){
	chat = Y / fL;
	fhat = f;
	Qhat = Q;
      }
    }
  }
  //now refine the estimators using Newton Raphson
  refine(y,p);
  
}

///Tolerance for Newton Raphson optimiser
#define NEWTON_TOL 1e-10
///Maximum number of Newton Raphson iterations
#define NEWTON_MAX_ITR 10
/// alpha=0: max log(periodogram), alpha=1: max standard periodogram
#define ALPHA 0

///Refine the Doppler estimate using Newton's method
///For details, see Andre Pollok's ASRP workbook, page 34-35.
void CoherentMackenthunWithDoppler::refine(const std::vector<complex>& y, const std::vector<complex>& p)
{
  
  hardDecisionsAndDerotate(y,p); //output goes to working memory variable yf
  
  MSGTYPE newtonstep = NEWTON_TOL + 1.0; //something bigger than NEWTONTOL to start with
  int numiters = 0; //number of iterations performed
  //run Newton's method
  complex Y(0,0);
  while( (std::fabs(newtonstep) > NEWTON_TOL) && (numiters < NEWTON_MAX_ITR) ){
    
    Y = complex(0,0);
    complex Ydf(0,0);
    complex Ydf2(0,0);
    for(unsigned int i = 0; i < L; i++) {
      int n = PUD[i];
      X[i] = yf[n] * std::polar<MSGTYPE>(1.0, -2*PI*fhat*Ts*(n+1)); //setup reused vector X
      Y += X[i]; ///compute Y (beware, this notation is Y/L in the estimate function!)
      Ydf += Ts * (n+1) * X[i]; ///compute first derivative of Y
      Ydf2 += Ts*Ts * (n+1)*(n+1) * X[i];
    }
    Y = complex(1.0/L,0.0) * Y; //normalise
    Ydf =  complex(0.0, -2*PI/L) * Ydf;
    Ydf2 = complex(-4.0*PI*PI/L,0.0) * Ydf2;

    MSGTYPE I = std::norm(Y); //compute I 
    MSGTYPE Idf = 2.0*std::real( Ydf * std::conj(Y) ); //first derivative of I
    MSGTYPE Idf2 = 2.0*std::real(  Ydf2 * std::conj(Y) ) + 2.0*std::norm(Ydf); //second derivative 
    MSGTYPE G = Idf*Idf;
    MSGTYPE H = Idf2;

    //Newton iterate with monotonic function kappa(I) = I^alpha.
    //For ALPHA=1, this is the standard Newton iterate inv(H)*[Idf; Idfr].
    newtonstep = Idf / ((ALPHA-1)/I*G + H);
    fhat = fhat - newtonstep; //update fhat    
    numiters += 1;
    
  }
  //get the update phase estimate
  chat = Y;

}

///Given the existing phase and Doppler estimates ahat and fhat, this computes hard decisions for the data and derotates the 
///recieved vector y in preparations for the refine method
void CoherentMackenthunWithDoppler::hardDecisionsAndDerotate(const std::vector<complex>& y, const std::vector<complex>& p)
{

  //data indices first
  MSGTYPE argchat = std::arg<MSGTYPE>(chat);
  for(unsigned int i = 0; i < D.size(); i++){
    int n = D[i];
    MSGTYPE u = round(M/(2*PI)*(std::arg<MSGTYPE>(y[n]) - argchat - 2*PI*(fhat*Ts*(n+1)))); // hard decision
    yf[n] = y[n] * std::polar<MSGTYPE>(1.0, -2*PI*u/M);
  }
  //now pilot indices
  for(unsigned int i = 0; i < P.size(); i++){
    int n = P[i];
    yf[n] = y[n] * std::conj<MSGTYPE>(p[i]);
  }

}


CoherentMackenthunWithDopplerAndDopplerRate::CoherentMackenthunWithDopplerAndDopplerRate(
							     const std::vector<int>& D,
							     const std::vector<int>& P,
							     const std::vector<complex>& p,
							     const unsigned int M,
							     const MSGTYPE fmin,
							     const MSGTYPE fmax,
							     const MSGTYPE frmin_,
							     const MSGTYPE frmax_,
							     const MSGTYPE T,
							     const MSGTYPE search_oversample) 
  : CoherentMackenthunWithDoppler(D,P,p,M,fmin,fmax,T,search_oversample),
    frmin(frmin_),
    frmax(frmax_)
{
  int indif = max_index-min_index;
  int frtrials = (int)ceil(search_oversample*M*T*T*indif*indif*(frmax - frmin));
  if( frtrials == 0 ) frtrials = 1; //atleast one trial
  frstep = (frmax - frmin)/frtrials;
  if( frstep == 0.0) frstep = 1.0; //prevent potential infinite loop
  //std::cout << "grid size = " << (frtrials * (int)ceil(search_oversample*M*T*indif*(fmax - fmin))) << std::endl;
}

/// Run the estimator
void CoherentMackenthunWithDopplerAndDopplerRate::estimate(const std::vector<complex>& y, const std::vector<complex>& p)
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
  
  MSGTYPE Qhat = -1.0;
  chat = complex(0,0);
  for(MSGTYPE f = fmin; f <= fmax; f += fstep){
    for(MSGTYPE fr = frmin; fr <= frmax; fr += frstep){
      
      //setup sequences and sort
      complex Y(0.0,0.0);
      for(unsigned int i=0; i < P.size(); i++) 
	Y += y[P[i]] * std::conj<MSGTYPE>(p[i]) * std::polar<MSGTYPE>(1.0, -2*PI*( f*Ts*(P[i]+1) + fr*Ts*Ts*(P[i]+1)*(P[i]+1) ) );
      for(unsigned int i=0; i < D.size(); i++){
	int n = z[i].i; //use last permutation order (should make sort faster)
	//int n = i;
	//Use D[n]+1 here and above to be consistent with MATLAB indexing
	complex DCtemp = std::polar<MSGTYPE>(1.0, -2*PI*( f*Ts*(D[n]+1) + fr*Ts*Ts*(D[n]+1)*(D[n]+1) ) ); 
	MSGTYPE phi = std::arg( y[D[n]] * DCtemp ); 
	MSGTYPE u = w*round(phi/w);
	g[n] = y[D[n]] * std::polar<MSGTYPE>(1,-u) * DCtemp;
	z[i] = IndexedReal(phi-u, n);
	Y += g[n];
      }
      sort(z.begin(), z.end());

      MSGTYPE fL = (MSGTYPE) L;
      for(unsigned int k = 0; k < M*D.size(); k++ ){
	Y += nu*g[sigma(k)];
	g[sigma(k)] *= eta;
	MSGTYPE Q = norm(Y) / fL;
	if( Q > Qhat ){
	  chat = Y / fL;
	  fhat = f;
	  frhat = fr;
	  Qhat = Q;
	}
      }
    }
  }
  //now refine the estimators using Newton Raphson
  refine(y,p);
  
}

void CoherentMackenthunWithDopplerAndDopplerRate::refine(const std::vector<complex>& y, const std::vector<complex>& p) {
  
  hardDecisionsAndDerotate(y,p); //output goes to working memory variable yf
  
  MSGTYPE newtonstep = NEWTON_TOL + 1.0; //something bigger than NEWTONTOL to start with
  int numiters = 0; //number of iterations performed
  //run Newton's method
  complex Y(0,0);
  
  while( (std::fabs(newtonstep) > NEWTON_TOL) && (numiters < NEWTON_MAX_ITR) ){
    
    //setup reused vector X
    Y = complex(0,0);
    complex Ydf(0,0);
    complex Ydfr(0,0);
    complex Ydf2(0,0);
    complex Ydfr2(0,0);
    complex Ydfrdf(0,0);
    for(unsigned int i = 0; i < L; i++) { 
      int n = PUD[i];
      X[i] = yf[n] * std::polar<MSGTYPE>(1.0, -2*PI*( fhat*Ts*(n+1) + frhat*Ts*Ts*(n+1)*(n+1) ));
      Y += X[i]; ///compute Y (beware, this notation is Y/L in the estimate function!)
      Ydf += Ts * (n+1) * X[i]; /// first derivative of Y with respect to f
      Ydfr += Ts*Ts*(n+1)*(n+1) * X[i]; /// first derivative of Y with respect to fr
      Ydf2 += Ts*Ts * (n+1)*(n+1) * X[i];  // second derivative of Y with respect to f
      Ydfr2 += std::pow(Ts*(n+1), (MSGTYPE)4.0) * X[i]; // second derivative of Y with respect to fr
      Ydfrdf += std::pow(Ts*(n+1), (MSGTYPE)3.0) * X[i]; // derivative of Y with respect to fr and f
    }
    Y = complex(1.0/L,0.0) * Y; //normalise    
    Ydf =  complex(0.0, -2*PI/L) * Ydf;
    Ydfr = complex(0.0, -2*PI/L) * Ydfr;
    Ydf2 = complex(-4.0*PI*PI/L,0.0) * Ydf2;
    Ydfr2 = complex(-4.0*PI*PI/L,0.0) * Ydfr2;
    Ydfrdf = complex(-4.0*PI*PI/L,0.0) * Ydfrdf;

    MSGTYPE I = std::norm(Y); //compute I 
    MSGTYPE Idf = 2.0*std::real( Ydf * std::conj(Y) ); //first derivative of I with respect to f
    MSGTYPE Idf2 = 2.0*std::real(  Ydf2 * std::conj(Y) ) + 2.0*std::norm(Ydf); //second derivative with respect to f
    MSGTYPE Idfr = 2.0*std::real( Ydfr * std::conj(Y) ); //first derivative with respect to fr
    MSGTYPE Idfr2 = 2.0*std::real(  Ydfr2 * std::conj(Y) ) + 2.0*std::norm(Ydfr); //second derivative with respect to fr
    MSGTYPE Idfrdf = 2.0*std::real(  Ydfrdf * std::conj(Y)  + Ydf * std::conj(Ydfr) ); //derivative with respect to fr and f
    
    //compute Hessian matrix and it's 2 by 2 inverse
    MSGTYPE M11 = (ALPHA-1.0)/I * Idf*Idf + Idf2;
    MSGTYPE M12 = (ALPHA-1.0)/I * Idf*Idfr + Idfrdf;
    MSGTYPE M21 = M12;
    MSGTYPE M22 = (ALPHA-1.0)/I * Idfr*Idfr + Idfr2;
    MSGTYPE det = M11*M22 - M12*M21; //determinant
    MSGTYPE I11 = M22/det;
    MSGTYPE I12 = -M12/det;
    MSGTYPE I21 = -M21/det;
    MSGTYPE I22 = M11/det;

    //Newton iterate with monotonic function kappa(I) = I^alpha.
    //For ALPHA=1, this is the standard Newton iterate inv(H)*[Idf; Idfr].
    MSGTYPE newtonstepf = I11*Idf + I12*Idfr;
    MSGTYPE newtonstepfr = I21*Idf + I22*Idfr;
    fhat = fhat - newtonstepf; //update fhat    
    frhat = frhat - newtonstepfr; //update frhat    
    numiters += 1;
    
    //step size for termination
    newtonstep = std::fabs(newtonstepf) + std::fabs(newtonstepfr);

  }
  //get the update phase estimate
  chat = Y;

}

void CoherentMackenthunWithDopplerAndDopplerRate::hardDecisionsAndDerotate(const std::vector<complex>& y, const std::vector<complex>& p){
  
  //data indices first
  MSGTYPE argchat = std::arg<MSGTYPE>(chat);
  for(unsigned int i = 0; i < D.size(); i++){
    int n = D[i];
    MSGTYPE u = round(M/(2*PI)*(std::arg<MSGTYPE>(y[n]) - argchat - 2*PI*( fhat*Ts*(n+1) + frhat*Ts*Ts*(n+1)*(n+1)))); // hard decision
    yf[n] = y[n] * std::polar<MSGTYPE>(1.0, -2*PI*u/M);
  }
  //now pilot indices
  for(unsigned int i = 0; i < P.size(); i++){
    int n = P[i];
    yf[n] = y[n] * std::conj<MSGTYPE>(p[i]);
  }
}
