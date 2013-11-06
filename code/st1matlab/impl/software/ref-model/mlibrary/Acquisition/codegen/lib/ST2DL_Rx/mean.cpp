/*
 * mean.cpp
 *
 * Code generation for function 'mean'
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
#include "mean.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 *
 */
creal_T mean(const creal_T x[3990])
{
  creal_T y;
  int32_T k;
  y = x[0];
  for (k = 0; k < 3989; k++) {
    y.re += x[k + 1].re;
    y.im += x[k + 1].im;
  }

  return eml_div(y, 3990.0);
}

/* End of code generation (mean.cpp) */
