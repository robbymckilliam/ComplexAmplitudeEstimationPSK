/*
 * sum.cpp
 *
 * Code generation for function 'sum'
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
#include "sum.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 *
 */
creal_T sum(const emxArray_creal_T *x)
{
  creal_T y;
  int32_T vlen;
  int32_T k;
  if (x->size[1] == 0) {
    y.re = 0.0;
    y.im = 0.0;
  } else {
    vlen = x->size[1];
    y = x->data[0];
    for (k = 2; k <= vlen; k++) {
      y.re += x->data[k - 1].re;
      y.im += x->data[k - 1].im;
    }
  }

  return y;
}

/* End of code generation (sum.cpp) */
