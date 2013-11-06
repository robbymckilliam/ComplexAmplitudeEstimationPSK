/*
 * dcblockST2DL_emxAPI.cpp
 *
 * Code generation for function 'dcblockST2DL_emxAPI'
 *
 * C source code generated on: Thu Nov  8 13:52:46 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dcblockST2DL.h"
#include "freqest_periodogram_coder.h"
#include "mf_coder.h"
#include "dcblockST2DL_emxAPI.h"
#include "dcblockST2DL_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
emxArray_creal_T *emxCreateND_creal_T(int32_T numDimensions, int32_T *size)
{
  emxArray_creal_T *emx;
  int32_T numEl;
  int32_T loop_ub;
  int32_T i;
  c_emxInit_creal_T(&emx, numDimensions);
  numEl = 1;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (creal_T *)calloc((uint32_T)numEl, sizeof(creal_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_creal_T *emxCreateWrapperND_creal_T(creal_T *data, int32_T
  numDimensions, int32_T *size)
{
  emxArray_creal_T *emx;
  int32_T numEl;
  int32_T loop_ub;
  int32_T i;
  c_emxInit_creal_T(&emx, numDimensions);
  numEl = 1;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_creal_T *emxCreateWrapper_creal_T(creal_T *data, int32_T rows, int32_T
  cols)
{
  emxArray_creal_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  c_emxInit_creal_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_creal_T *emxCreate_creal_T(int32_T rows, int32_T cols)
{
  emxArray_creal_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  c_emxInit_creal_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (creal_T *)calloc((uint32_T)numEl, sizeof(creal_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

void emxDestroyArray_creal_T(emxArray_creal_T *emxArray)
{
  emxFree_creal_T(&emxArray);
}

/* End of code generation (dcblockST2DL_emxAPI.cpp) */
