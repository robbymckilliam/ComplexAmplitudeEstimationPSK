/*
 * abs.cpp
 *
 * Code generation for function 'abs'
 *
 * C source code generated on: Thu Nov  8 13:52:45 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dcblockST2DL.h"
#include "freqest_periodogram_coder.h"
#include "mf_coder.h"
#include "abs.h"
#include "dcblockST2DL_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  uint32_T uv2[2];
  int32_T k;
  for (k = 0; k < 2; k++) {
    uv2[k] = (uint32_T)x->size[k];
  }

  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int32_T)uv2[1];
  emxEnsureCapacity((emxArray__common *)y, k, (int32_T)sizeof(real_T));
  for (k = 0; k <= x->size[1] - 1; k++) {
    y->data[k] = fabs(x->data[k]);
  }
}

/* End of code generation (abs.cpp) */
