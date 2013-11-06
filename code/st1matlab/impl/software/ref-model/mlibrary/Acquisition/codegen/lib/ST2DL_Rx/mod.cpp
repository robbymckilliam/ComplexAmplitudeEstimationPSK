/*
 * mod.cpp
 *
 * Code generation for function 'mod'
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
#include "mod.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 *
 */
void b_mod(const real_T x[3960], real_T y, real_T r[3960])
{
  int32_T k;
  real_T b_r;
  real_T b_x;
  for (k = 0; k < 3960; k++) {
    if (y == 0.0) {
      b_r = x[k];
    } else if (y == floor(y)) {
      b_r = x[k] - floor(x[k] / y) * y;
    } else {
      b_r = x[k] / y;
      if (fabs(b_r) > 4.503599627370496E+15) {
        b_x = b_r;
      } else if (b_r >= 0.5) {
        b_x = floor(b_r + 0.5);
      } else if (b_r > -0.5) {
        b_x = -0.0;
      } else {
        b_x = ceil(b_r - 0.5);
      }

      if (fabs(b_r - b_x) <= 2.2204460492503131E-16 * fabs(b_r)) {
        b_r = 0.0;
      } else {
        b_r = (b_r - floor(b_r)) * y;
      }
    }

    r[k] = b_r;
  }
}

/* End of code generation (mod.cpp) */
