/*
 * power.cpp
 *
 * Code generation for function 'power'
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
#include "power.h"
#include "ST2DL_Rx_emxutil.h"
#include "ST2DL_Rx_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 *
 */
void power(const emxArray_real_T *a, emxArray_real_T *y)
{
  uint32_T uv1[2];
  int32_T i6;
  int32_T k;
  for (i6 = 0; i6 < 2; i6++) {
    uv1[i6] = (uint32_T)a->size[i6];
  }

  i6 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int32_T)uv1[1];
  emxEnsureCapacity((emxArray__common *)y, i6, (int32_T)sizeof(real_T));
  i6 = y->size[1];
  for (k = 0; k <= i6 - 1; k++) {
    y->data[k] = rt_powd_snf(a->data[k], 2.0);
  }
}

/* End of code generation (power.cpp) */
