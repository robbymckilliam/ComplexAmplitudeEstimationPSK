/*
 * dcblockST2DL_emxAPI.h
 *
 * Code generation for function 'dcblockST2DL_emxAPI'
 *
 * C source code generated on: Thu Nov  8 13:52:46 2012
 *
 */

#ifndef __DCBLOCKST2DL_EMXAPI_H__
#define __DCBLOCKST2DL_EMXAPI_H__
/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "dcblockST2DL_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern emxArray_creal_T *emxCreateND_creal_T(int32_T numDimensions, int32_T *size);
extern emxArray_creal_T *emxCreateWrapperND_creal_T(creal_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_creal_T *emxCreateWrapper_creal_T(creal_T *data, int32_T rows, int32_T cols);
extern emxArray_creal_T *emxCreate_creal_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_creal_T(emxArray_creal_T *emxArray);
#endif
/* End of code generation (dcblockST2DL_emxAPI.h) */
