/*
 * round.cpp
 *
 * Code generation for function 'round'
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
#include "round.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 *
 */
void b_round(real_T x[3960])
{
  int32_T k;
  real_T b_x;
  for (k = 0; k < 3960; k++) {
    if (fabs(x[k]) > 4.503599627370496E+15) {
      b_x = x[k];
    } else if (x[k] >= 0.5) {
      b_x = floor(x[k] + 0.5);
    } else if (x[k] > -0.5) {
      b_x = x[k] * 0.0;
    } else {
      b_x = ceil(x[k] - 0.5);
    }

    x[k] = b_x;
  }
}

/* End of code generation (round.cpp) */
