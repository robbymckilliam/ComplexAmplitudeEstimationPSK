/* 
 * File:   MainLoop.h
 * Author: mckillrg
 *
 * Created on 4 March 2013, 3:13 PM
 */

#ifndef MAINLOOP_H
#define	MAINLOOP_H

#include <complex>

#ifndef VALTYPE
#define VALTYPE double
#endif
#ifndef MSGTYPE
#define MSGTYPE VALTYPE
#endif

using namespace std;

namespace TimeOffset {

#ifdef __ARM_NEON__
  #warning "Compiling with NEON instruction set"
  #include "arm_neon.h"
  #include "neon_mathfun.h"
    //some useful macros
#define RMEMSTACKED(m) {std::real(rmem[m]), std::imag(rmem[m]), std::real(rmem[m+mstep]), std::imag(rmem[m+mstep])};
#define RMEMFILL(v,m) v[0]=std::real(rmem[m]); v[1]=std::imag(rmem[m]); v[2]=std::real(rmem[m+mstep]); v[3]=std::imag(rmem[m+mstep]);
#define PULSESTACKED(i) {pulsetable[i], pulsetable[i], pulsetable[i+istep], pulsetable[i+istep]};
#define PULSEFILL(v,i) v[0]=pulsetable[i]; v[1]=pulsetable[i]; v[2]=pulsetable[i+istep]; v[3]=pulsetable[i+istep];
     /** Runs the main filter loop, this is where most time is spent */
    inline complex<VALTYPE> mainFilterLoop(int mfrom, int mto, int mstep, int ifrom, int istep, const complex<VALTYPE>* rmem, const VALTYPE* pulsetable) {

      float32x4_t sum;
      int numsteps = (mto - mfrom + 1)/mstep;
      if( numsteps%2 == 0) { //load zeros if even
	float32_t sump[4] = {0,0,0,0};
	sum = vld1q_f32(sump);
      } else { //load first one and increment counters if it's odd
	complex<VALTYPE> first = rmem[mfrom] * pulsetable[ifrom];
	float32_t sump[4] = {std::real(first),std::imag(first),0,0};
	sum = vld1q_f32(sump);
	mfrom += mstep;
	ifrom += istep;
      }
      int istep2 = istep*2, mstep2 = mstep*2, mstep8 = mstep*8, istep8 = istep*8; //new counters
      
      float32_t ap[4], bp[4];
      //main loop
      #pragma unroll 
      while(mfrom <= mto) {
	RMEMFILL(ap,mfrom);
	PULSEFILL(bp,ifrom);
	float32x4_t a = vld1q_f32(ap);
	float32x4_t b = vld1q_f32(bp);
	sum = vmlaq_f32(sum,a,b); // 4x multiply accumulate
	mfrom += mstep2;
	ifrom += istep2;
      }
      
      float32_t sump[4];
      vst1q_f32(sump, sum);
      return complex<VALTYPE>(sump[0] + sump[2],sump[1] + sump[3]);
    }
#else
    /** Runs the main filter loop, this is where most time is spent */
    static inline complex<VALTYPE> mainFilterLoop(int mfrom, int mto, int mstep, int ifrom, int istep, const complex<VALTYPE>* rmem, const  VALTYPE* pulsetable) {
        complex<VALTYPE> sum(0, 0);
        
        //mainloop, this will greatly benefit from a mac instruction
        for(int m = mfrom, i = ifrom; m<=mto; m+=mstep, i+=istep) 
            sum += rmem[m] * pulsetable[i];   
        
        return sum;
    }
#endif

   /** 
    * Runs the main filter loop. This assumes you have arranged memory so that the pulsetable 
    * can be stepped through with unit increments.
    */
   static inline complex<VALTYPE> mainBankedFilterLoop(int mfrom, int mto, int mstep, const complex<VALTYPE>* rmem, const VALTYPE* pulsetable) {
       return mainFilterLoop(mfrom,mto,mstep,0,1,rmem,pulsetable);
   }
    
}

#endif	/* MAINLOOP_H */

