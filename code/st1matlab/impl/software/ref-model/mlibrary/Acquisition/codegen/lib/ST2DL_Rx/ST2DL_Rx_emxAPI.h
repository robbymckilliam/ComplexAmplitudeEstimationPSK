/*
 * ST2DL_Rx_emxAPI.h
 *
 * Code generation for function 'ST2DL_Rx_emxAPI'
 *
 * C source code generated on: Fri Nov 23 11:11:00 2012
 *
 */

#ifndef __ST2DL_RX_EMXAPI_H__
#define __ST2DL_RX_EMXAPI_H__
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
extern emxArray_creal_T *emxCreateND_creal_T(int32_T numDimensions, int32_T *size);
extern emxArray_creal_T *emxCreateWrapperND_creal_T(creal_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_creal_T *emxCreateWrapper_creal_T(creal_T *data, int32_T rows, int32_T cols);
extern emxArray_creal_T *emxCreate_creal_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_creal_T(emxArray_creal_T *emxArray);
#endif
/* End of code generation (ST2DL_Rx_emxAPI.h) */
