/*
 * freqest_periodogram_coder.h
 *
 * Code generation for function 'freqest_periodogram_coder'
 *
 * C source code generated on: Fri Nov 23 11:10:59 2012
 *
 */

#ifndef __FREQEST_PERIODOGRAM_CODER_H__
#define __FREQEST_PERIODOGRAM_CODER_H__
/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "ST2DL_Rx_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern creal_T eml_div(const creal_T x, real_T y);
extern void freqest_periodogram_coder(const creal_T RxSamp[16000], const real_T tSamp[16000], real_T T, real_T OS, real_T Lsin, real_T taumin, real_T taumax, real_T *fhat, real_T *phihat);
#endif
/* End of code generation (freqest_periodogram_coder.h) */
