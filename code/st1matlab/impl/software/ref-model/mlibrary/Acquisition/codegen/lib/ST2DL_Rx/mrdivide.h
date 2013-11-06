/*
 * mrdivide.h
 *
 * Code generation for function 'mrdivide'
 *
 * C source code generated on: Fri Nov 23 11:10:59 2012
 *
 */

#ifndef __MRDIVIDE_H__
#define __MRDIVIDE_H__
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
extern void b_mrdivide(const real_T A[3960], real_T B, real_T y[3960]);
extern creal_T c_mrdivide(const creal_T A, const creal_T B);
extern real_T mrdivide(real_T A, real_T B);
#endif
/* End of code generation (mrdivide.h) */
