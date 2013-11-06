/*
 * ipermute.cpp
 *
 * Code generation for function 'ipermute'
 *
 * C source code generated on: Fri Nov 23 11:10:59 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "EsNoest_coder.h"
#include "dcblockST2DL.h"
#include "freqest_periodogram_coder.h"
#include "mf_coder.h"
#include "ipermute.h"
#include "ST2DL_Rx_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 *
 */
void ipermute(const emxArray_creal_T *b, emxArray_creal_T *a)
{
  int32_T i5;
  int32_T unnamed_idx_1;
  i5 = a->size[0] * a->size[1];
  a->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)a, i5, (int32_T)sizeof(creal_T));
  unnamed_idx_1 = b->size[0];
  i5 = a->size[0] * a->size[1];
  a->size[1] = unnamed_idx_1;
  emxEnsureCapacity((emxArray__common *)a, i5, (int32_T)sizeof(creal_T));
  unnamed_idx_1 = b->size[0] - 1;
  for (i5 = 0; i5 <= unnamed_idx_1; i5++) {
    a->data[i5] = b->data[i5];
  }
}

/* End of code generation (ipermute.cpp) */
