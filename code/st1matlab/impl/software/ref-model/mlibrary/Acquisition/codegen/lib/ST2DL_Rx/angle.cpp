/*
 * angle.cpp
 *
 * Code generation for function 'angle'
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
#include "angle.h"
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
void angle(const creal_T x[3960], real_T y[3960])
{
  int32_T k;
  for (k = 0; k < 3960; k++) {
    y[k] = rt_atan2d_snf(x[k].im, x[k].re);
  }
}

/* End of code generation (angle.cpp) */
