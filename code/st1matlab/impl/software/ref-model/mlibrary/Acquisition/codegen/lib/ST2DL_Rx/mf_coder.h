/*
 * mf_coder.h
 *
 * Code generation for function 'mf_coder'
 *
 * C source code generated on: Fri Nov 23 11:10:59 2012
 *
 */

#ifndef __MF_CODER_H__
#define __MF_CODER_H__
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
extern void mf_coder(const creal_T rx[16000], real_T tauhat, real_T rolloff, real_T q, real_T T, const real_T P[30], const real_T D[3960], real_T RrcNsym, emxArray_creal_T *mfout);
#endif
/* End of code generation (mf_coder.h) */
