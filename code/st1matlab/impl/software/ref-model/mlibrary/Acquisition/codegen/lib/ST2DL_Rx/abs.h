/*
 * abs.h
 *
 * Code generation for function 'abs'
 *
 * C source code generated on: Fri Nov 23 11:10:59 2012
 *
 */

#ifndef __ABS_H__
#define __ABS_H__
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
extern real_T b_abs(const creal_T x);
extern void c_abs(const emxArray_real_T *x, emxArray_real_T *y);
#endif
/* End of code generation (abs.h) */
