/*
 * rdivide.cpp
 *
 * Code generation for function 'rdivide'
 *
 * C source code generated on: Thu Nov  8 13:52:45 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dcblockST2DL.h"
#include "freqest_periodogram_coder.h"
#include "mf_coder.h"
#include "rdivide.h"
#include "dcblockST2DL_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void rdivide(const emxArray_real_T *x, const emxArray_real_T *y, emxArray_real_T
             *z)
{
  int32_T i8;
  int32_T loop_ub;
  i8 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)z, i8, (int32_T)sizeof(real_T));
  loop_ub = x->size[0] * x->size[1] - 1;
  for (i8 = 0; i8 <= loop_ub; i8++) {
    z->data[i8] = x->data[i8] / y->data[i8];
  }
}

/* End of code generation (rdivide.cpp) */
