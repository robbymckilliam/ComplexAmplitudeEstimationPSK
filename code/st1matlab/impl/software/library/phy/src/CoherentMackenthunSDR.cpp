#include "CoherentMackenthunSDR.h"
#include "utilities.h"
//#include "mex.h"
#include <math.h>
#include <float.h>
#include <algorithm>
#include <sstream> 
#include <iostream>

using namespace std;
 
//#define PI 	3.1415926535897932384626433
#define PI 	3.1415926f

#ifdef _WIN32
#define log2(x) 			(log(x)/log((CMACK_MSGTYPE)2.0))
#define round(x)			(floor((x)+0.5))
#endif

/// \brief	Standard constructor.
CoherentMackenthun::CoherentMackenthun(
 				       const std::vector<int>* Din, 
 				       const std::vector<int>* Pin,
 				       const std::vector<mycomplex>* pin,
  				       const int Min) :
  L(Din->size() + Pin->size()),
  M(Min),
  Dptr(Din),
  Pptr(Pin),
  pptr(pin),
  absD(Din->size()),
  absP(Pin->size()),
  w(2*PI/M),
  //eta(std::polar<CMACK_MSGTYPE>(1.0,w)),
  //nu(std::polar<CMACK_MSGTYPE>(1.0,w) - mycomplex(1.0,0.0))
  eta(mycomplex(0.0,1.0)),
  nu(mycomplex(0.0,1.0) - mycomplex(1.0,0.0))
{
    //cout << eta << " " << std::polar<CMACK_MSGTYPE>(1.0,w) << endl;
    //cout << nu << " " << std::polar<CMACK_MSGTYPE>(1.0,w) - mycomplex(1.0,0.0) << endl;
  if( pin->size() != Pin->size() ) throw std::string("The number of pilots symbols does not match the number of indices in Pin");
  if( M < 1 )  throw std::string("M must be a positive integer");
  //allocate memory
  z.resize(Din->size(), IndexedReal(0,0));
  for(int i = 0; i < z.size(); i++) z[i] = IndexedReal(0,i); //initialise z with the identity permutation 
  g.resize(Din->size(), mycomplex(0,0)); 
}

/// \brief	Standard destructor.
///
CoherentMackenthun::~CoherentMackenthun(void)
{
  delete Dptr;
  delete Pptr;
  delete pptr;
  z.clear();
  g.clear();
}

/// Run the estimator
/// This version will delete the data vector when finished. If you don't want this, then use the & version. 
void CoherentMackenthun::estimate(const std::vector<mycomplex>* y)
{
  estimate(*y);
  delete y;
}

/// Run the estimator
/// This version will delete the data vector when finished. If you don't want this, then use the & version. 
void CoherentMackenthun::estimate(const std::vector<mycomplex>* y, const std::vector<mycomplex>* p)
{
  estimate(*y, *p);
  delete y;
}

/// Run the estimator on data y
void CoherentMackenthun::estimate(const std::vector<mycomplex>& y)
{
  estimate(y, *pptr);
}

/// Run the estimator with data y and pilot symbols p
void CoherentMackenthun::estimate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p)
{
  if( y.size() != L ) {
    std::ostringstream str;
    str << "Exception in function CoherentMackenthun::estimate. The data length y is ";
    str << y.size() << " but the specified symbol length L is " << L;
    throw str.str();
  }
  if( p.size() != Pptr->size() ) {
    std::ostringstream str;
    str << "Exception in function CoherentMackenthun::estimate. The number of pilot symbols given is ";
    str << p.size() << " but the number of indices is " << Pptr->size();
    throw str.str();
  }

  //references for convenience
  const std::vector<int>& P = *Pptr;
  const std::vector<int>& D = *Dptr;

  CMACK_MSGTYPE w_inv = M/(2*PI);
  //setup sequences and sort
  mycomplex Y(0.0,0.0);
  for(int i=0; i < P.size(); i++) Y += y[P[i]]*conj(p[i]);
  for(int i=0; i < D.size(); i++){
    CMACK_MSGTYPE phi = std::arg(y[D[i]]);
    //CMACK_MSGTYPE u = w*round(phi/w);
    CMACK_MSGTYPE u = w*round(phi*w_inv);
    g[i] = y[D[i]] * std::polar<CMACK_MSGTYPE>(1,-u);
    z[i] = IndexedReal(phi-u, i);
    Y += g[i];
  }
 
  sort(z.begin(), z.end());
  /*
  CMACK_MSGTYPE fL = (CMACK_MSGTYPE) L;
  chat = Y / fL;
  Qhat = norm(Y) / fL;
  */
  CMACK_MSGTYPE fL = 1/(CMACK_MSGTYPE) L;
  chat = Y * fL;
  Qhat = norm(Y) * fL;
  int sigma_k;
      
  for( int k = 0; k < M*D.size(); k++ ){
    //Y += nu*g[sigma(k)];
    //g[sigma(k)] *= eta;
    //CMACK_MSGTYPE Q = norm(Y) / fL;
      sigma_k = sigma(k);
      Y += nu*g[sigma_k];
      g[sigma_k] *= eta;
      CMACK_MSGTYPE Q = norm(Y) * fL;
    if( Q > Qhat ){
	chat = Y * fL;
	Qhat = Q;
    }
   }
}

//SEARCH_OVERSAMPLE adjusts how fine the search grid is (smaller is faster, larger is more accurate)
//This probably shouldn't be hardcoded like this!
#define SEARCH_OVERSAMPLE 2.0

CoherentMackenthunWithDoppler::CoherentMackenthunWithDoppler(
							     const std::vector<int>* D,
							     const std::vector<int>* P,
							     const std::vector<mycomplex>* p,
							     const int M,
							     const CMACK_MSGTYPE fmin_,
							     const CMACK_MSGTYPE fmax_,
							     const CMACK_MSGTYPE T) 
  : CoherentMackenthun(D,P,p,M), 
    fmin(fmin_),
    fmax(fmax_),
    Ts(T)
{
  int ftrials = (int)ceil(SEARCH_OVERSAMPLE*M*(T*L)*(fmax - fmin));
  if( ftrials == 0 ) ftrials = 1; //atleast one trial
  fstep = (fmax - fmin)/ftrials;
  if( fstep == 0.0) fstep = 1.0; //prevent potential infinite loop
  //allocate working memory
  yf.resize(L, mycomplex(0,0));
  X.resize(L, mycomplex(0,0));
}

/// Run the estimator
void CoherentMackenthunWithDoppler::estimate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p)
{
  if( y.size() != L ) {
    std::ostringstream str;
    str << "Exception in function CoherentMackenthun::estimate. The data length y is ";
    str << y.size() << " but the specified symbol length L is " << L;
    throw str.str();
  }
  if( p.size() != Pptr->size() ) {
    std::ostringstream str;
    str << "Exception in function CoherentMackenthun::estimate. The number of pilot symbols given is ";
    str << p.size() << " but the number of indices is " << Pptr->size();
    throw str.str();
 }

  //references for convenience
  const std::vector<int>& P = *Pptr;
  const std::vector<int>& D = *Dptr;
  
  CMACK_MSGTYPE Qhat = -1.0;
  chat = mycomplex(0,0);
  for(CMACK_MSGTYPE f = fmin; f <= fmax; f += fstep){

    //setup sequences and sort
    mycomplex Y(0.0,0.0);
    for(int i=0; i < P.size(); i++) Y += y[P[i]] * conj(p[i]) * std::polar<CMACK_MSGTYPE>(1.0, -2*PI*f*Ts*(P[i]+1));
    for(int i=0; i < D.size(); i++){
      int n = z[i].i; //use last permutation order (should make sort faster)
      //int n = i;
      mycomplex DCtemp = std::polar<CMACK_MSGTYPE>(1.0, -2*PI*f*Ts*(D[n]+1));
      CMACK_MSGTYPE phi = std::arg( y[D[n]] * DCtemp ); //Use D[i]+1 here and above to be consistent with MATLAB indexing
      CMACK_MSGTYPE u = w*round(phi/w);
      g[n] = y[D[n]] * std::polar<CMACK_MSGTYPE>(1,-u) * DCtemp;
      z[i] = IndexedReal(phi-u, n);
      Y += g[n];
    }
    sort(z.begin(), z.end());

    CMACK_MSGTYPE fL = (CMACK_MSGTYPE) L;
    for( int k = 0; k < (M+1)*D.size(); k++ ){
      Y += nu*g[sigma(k)];
      g[sigma(k)] *= eta;
      CMACK_MSGTYPE Q = norm(Y) / fL;
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
void CoherentMackenthunWithDoppler::refine(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p)
{
  
  hardDecisionsAndDerotate(y,p); //output goes to working memory variable yf
  
  CMACK_MSGTYPE newtonstep = NEWTON_TOL + 1.0; //something bigger than NEWTONTOL to start with
  int numiters = 0; //number of iterations performed
  //run Newton's method
  mycomplex Y(0,0);
  while( (std::fabs(newtonstep) > NEWTON_TOL) && (numiters < NEWTON_MAX_ITR) ){
    
    //setup reused vector X
    for(int i = 0; i < L; i++) X[i] = yf[i] * std::polar<CMACK_MSGTYPE>(1.0, -2*PI*fhat*Ts*(i+1));

    ///compute Y (beware, this notation is Y/L in the estimate function!)
    Y = mycomplex(0,0);
    for(int i = 0; i < L; i++) Y += X[i];
    Y = mycomplex(1.0/L,0.0) * Y;
    
    ///compute first derivative of Y
    mycomplex Ydf(0,0);
    for(int i = 0; i < L; i++) Ydf += Ts * (i+1) * X[i];
    Ydf =  mycomplex(0.0, -2*PI/L) * Ydf;
    
    //compute second derivative of Y
    mycomplex Ydf2(0,0);
    for(int i = 0; i < L; i++) Ydf2 += Ts*Ts * (i+1)*(i+1) * X[i];
    Ydf2 = mycomplex(-4.0*PI*PI/L,0.0) * Ydf2;

    CMACK_MSGTYPE I = std::norm(Y); //compute I 
    CMACK_MSGTYPE Idf = 2.0*std::real( Ydf * std::conj(Y) ); //first derivative of I
    CMACK_MSGTYPE Idf2 = 2.0*std::real(  Ydf2 * std::conj(Y) ) + 2.0*std::norm(Ydf); //second derivative
    
    CMACK_MSGTYPE G = Idf*Idf;
    CMACK_MSGTYPE H = Idf2;

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
void CoherentMackenthunWithDoppler::hardDecisionsAndDerotate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p)
{
  //references for convenience
  const std::vector<int>& P = *Pptr;
  const std::vector<int>& D = *Dptr;

  //data indices first
  CMACK_MSGTYPE argchat = std::arg<CMACK_MSGTYPE>(chat);
  for(int i = 0; i < D.size(); i++){
    int n = D[i];
    CMACK_MSGTYPE u = round(M/(2*PI)*(std::arg<CMACK_MSGTYPE>(y[n]) - argchat - 2*PI*(fhat*Ts*(n+1)))); // hard decision
    yf[n] = y[n] * std::polar<CMACK_MSGTYPE>(1.0, -2*PI*u/M);
  }
  //now pilot indices
  for(int i = 0; i < P.size(); i++){
    int n = P[i];
    yf[n] = y[n] * std::conj<CMACK_MSGTYPE>(p[i]);
  }

}


CoherentMackenthunWithDopplerAndDopplerRate::CoherentMackenthunWithDopplerAndDopplerRate(
							     const std::vector<int>* D,
							     const std::vector<int>* P,
							     const std::vector<mycomplex>* p,
							     const int M,
							     const CMACK_MSGTYPE fmin,
							     const CMACK_MSGTYPE fmax,
							     const CMACK_MSGTYPE frmin_,
							     const CMACK_MSGTYPE frmax_,
							     const CMACK_MSGTYPE T) 
  : CoherentMackenthunWithDoppler(D,P,p,M,fmin,fmax,T),
    frmin(frmin_),
    frmax(frmax_)
{
  int frtrials = (int)ceil(SEARCH_OVERSAMPLE*M*T*T*L*L*(frmax - frmin));
  if( frtrials == 0 ) frtrials = 1; //atleast one trial
  frstep = (frmax - frmin)/frtrials;
  if( frstep == 0.0) frstep = 1.0; //prevent potential infinite loop
  
}

/// Run the estimator
void CoherentMackenthunWithDopplerAndDopplerRate::estimate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p)
{
  if( y.size() != L ) {
    std::ostringstream str;
    str << "Exception in function CoherentMackenthun::estimate. The data length y is ";
    str << y.size() << " but the specified symbol length L is " << L;
    throw str.str();
  }
  if( p.size() != Pptr->size() ) {
    std::ostringstream str;
    str << "Exception in function CoherentMackenthun::estimate. The number of pilot symbols given is ";
    str << p.size() << " but the number of indices is " << Pptr->size();
    throw str.str();
 }

  //references for convenience
  const std::vector<int>& P = *Pptr;
  const std::vector<int>& D = *Dptr;
  
  CMACK_MSGTYPE Qhat = -1.0;
  chat = mycomplex(0,0);
  for(CMACK_MSGTYPE f = fmin; f <= fmax; f += fstep){
    for(CMACK_MSGTYPE fr = frmin; fr <= frmax; fr += frstep){
      
      //setup sequences and sort
      mycomplex Y(0.0,0.0);
      for(int i=0; i < P.size(); i++) Y += y[P[i]] * std::conj<CMACK_MSGTYPE>(p[i]) * std::polar<CMACK_MSGTYPE>(1.0, -2*PI*( f*Ts*(P[i]+1) + fr*Ts*Ts*(P[i]+1)*(P[i]+1) ) );
      for(int i=0; i < D.size(); i++){
	int n = z[i].i; //use last permutation order (should make sort faster)
	//int n = i;
	//Use D[n]+1 here and above to be consistent with MATLAB indexing
	mycomplex DCtemp = std::polar<CMACK_MSGTYPE>(1.0, -2*PI*( f*Ts*(D[n]+1) + fr*Ts*Ts*(D[n]+1)*(D[n]+1) ) ); 
	CMACK_MSGTYPE phi = std::arg( y[D[n]] * DCtemp ); 
	CMACK_MSGTYPE u = w*round(phi/w);
	g[n] = y[D[n]] * std::polar<CMACK_MSGTYPE>(1,-u) * DCtemp;
	z[i] = IndexedReal(phi-u, n);
	Y += g[n];
      }
      sort(z.begin(), z.end());

      CMACK_MSGTYPE fL = (CMACK_MSGTYPE) L;
      for( int k = 0; k < (M+1)*D.size(); k++ ){
	Y += nu*g[sigma(k)];
	g[sigma(k)] *= eta;
	CMACK_MSGTYPE Q = norm(Y) / fL;
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

void CoherentMackenthunWithDopplerAndDopplerRate::refine(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p) {
  
  hardDecisionsAndDerotate(y,p); //output goes to working memory variable yf
  
  CMACK_MSGTYPE newtonstep = NEWTON_TOL + 1.0; //something bigger than NEWTONTOL to start with
  int numiters = 0; //number of iterations performed
  //run Newton's method
  mycomplex Y(0,0);
  
  while( (std::fabs(newtonstep) > NEWTON_TOL) && (numiters < NEWTON_MAX_ITR) ){
    
    //setup reused vector X
    for(int i = 0; i < L; i++) X[i] = yf[i] * std::polar<CMACK_MSGTYPE>(1.0, -2*PI*( fhat*Ts*(i+1) + frhat*Ts*Ts*(i+1)*(i+1)));

    ///compute Y (beware, this notation is Y/L in the estimate function!)
    Y = mycomplex(0,0);
    for(int i = 0; i < L; i++) Y += X[i];
    Y = mycomplex(1.0/L,0.0) * Y;
    
    ///compute first derivative of Y with respect to f
    mycomplex Ydf(0,0);
    for(int i = 0; i < L; i++) Ydf += Ts * (i+1) * X[i];
    Ydf =  mycomplex(0.0, -2*PI/L) * Ydf;

    ///compute first derivative of Y with respect to fr
    mycomplex Ydfr(0,0);
    for(int i = 0; i < L; i++) Ydfr += Ts*Ts*(i+1)*(i+1) * X[i];
    Ydfr = mycomplex(0.0, -2*PI/L) * Ydfr;
    
    //compute second derivative of Y with respect to f
    mycomplex Ydf2(0,0);
    for(int i = 0; i < L; i++) Ydf2 += Ts*Ts * (i+1)*(i+1) * X[i];
    Ydf2 = mycomplex(-4.0*PI*PI/L,0.0) * Ydf2;

    //compute second derivative of Y with respect to fr
    mycomplex Ydfr2(0,0);
    for(int i = 0; i < L; i++) Ydfr2 += std::pow(Ts*(i+1), (CMACK_MSGTYPE)4.0) * X[i];
    Ydfr2 = mycomplex(-4.0*PI*PI/L,0.0) * Ydfr2;

    //compute derivative of Y with respect to fr and f
    mycomplex Ydfrdf(0,0);
    for(int i = 0; i < L; i++) Ydfrdf += std::pow(Ts*(i+1), (CMACK_MSGTYPE)3.0) * X[i];
    Ydfrdf = mycomplex(-4.0*PI*PI/L,0.0) * Ydfrdf;

    CMACK_MSGTYPE I = std::norm(Y); //compute I 
    CMACK_MSGTYPE Idf = 2.0*std::real( Ydf * std::conj(Y) ); //first derivative of I with respect to f
    CMACK_MSGTYPE Idf2 = 2.0*std::real(  Ydf2 * std::conj(Y) ) + 2.0*std::norm(Ydf); //second derivative with respect to f
    CMACK_MSGTYPE Idfr = 2.0*std::real( Ydfr * std::conj(Y) ); //first derivative with respect to fr
    CMACK_MSGTYPE Idfr2 = 2.0*std::real(  Ydfr2 * std::conj(Y) ) + 2.0*std::norm(Ydfr); //second derivative with respect to fr
    CMACK_MSGTYPE Idfrdf = 2.0*std::real(  Ydfrdf * std::conj(Y)  + Ydf * std::conj(Ydfr) ); //derivative with respect to fr and f
    
    //compute Hessian matrix and it's 2 by 2 inverse
    CMACK_MSGTYPE M11 = (ALPHA-1.0)/I * Idf*Idf + Idf2;
    CMACK_MSGTYPE M12 = (ALPHA-1.0)/I * Idf*Idfr + Idfrdf;
    CMACK_MSGTYPE M21 = M12;
    CMACK_MSGTYPE M22 = (ALPHA-1.0)/I * Idfr*Idfr + Idfr2;
    CMACK_MSGTYPE det = M11*M22 - M12*M21; //determinant
    CMACK_MSGTYPE I11 = M22/det;
    CMACK_MSGTYPE I12 = -M12/det;
    CMACK_MSGTYPE I21 = -M21/det;
    CMACK_MSGTYPE I22 = M11/det;

    //Newton iterate with monotonic function kappa(I) = I^alpha.
    //For ALPHA=1, this is the standard Newton iterate inv(H)*[Idf; Idfr].
    CMACK_MSGTYPE newtonstepf = I11*Idf + I12*Idfr;
    CMACK_MSGTYPE newtonstepfr = I21*Idf + I22*Idfr;
    fhat = fhat - newtonstepf; //update fhat    
    frhat = frhat - newtonstepfr; //update frhat    
    numiters += 1;
    
    //step size for termination
    newtonstep = std::fabs(newtonstepf) + std::fabs(newtonstepfr);

  }
  //get the update phase estimate
  chat = Y;

}

void CoherentMackenthunWithDopplerAndDopplerRate::hardDecisionsAndDerotate(const std::vector<mycomplex>& y, const std::vector<mycomplex>& p){
  //references for convenience
  const std::vector<int>& P = *Pptr;
  const std::vector<int>& D = *Dptr;

  //data indices first
  CMACK_MSGTYPE argchat = std::arg<CMACK_MSGTYPE>(chat);
  for(int i = 0; i < D.size(); i++){
    int n = D[i];
    CMACK_MSGTYPE u = round(M/(2*PI)*(std::arg<CMACK_MSGTYPE>(y[n]) - argchat - 2*PI*( fhat*Ts*(n+1) + frhat*Ts*Ts*(n+1)*(n+1)))); // hard decision
    yf[n] = y[n] * std::polar<CMACK_MSGTYPE>(1.0, -2*PI*u/M);
  }
  //now pilot indices
  for(int i = 0; i < P.size(); i++){
    int n = P[i];
    yf[n] = y[n] * std::conj<CMACK_MSGTYPE>(p[i]);
  }
}
