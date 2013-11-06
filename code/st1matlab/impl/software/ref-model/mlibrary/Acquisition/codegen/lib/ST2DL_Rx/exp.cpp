/*
 * exp.cpp
 *
 * Code generation for function 'exp'
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
#include "exp.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void c_exp(creal_T x[3960]);

/* Function Definitions */

/*
 *
 */
static void c_exp(creal_T x[3960])
{
  int32_T k;
  real_T r;
  real_T x_im;
  real_T b_x_im;
  for (k = 0; k < 3960; k++) {
    r = exp(x[k].re / 2.0);
    x_im = x[k].im;
    b_x_im = x[k].im;
    x[k].re = r * (r * cos(x_im));
    x[k].im = r * (r * sin(b_x_im));
  }
}

/*
 *
 */
void b_exp(const creal_T x[3960], creal_T b_x[3960])
{
  memcpy((void *)&b_x[0], (void *)&x[0], 3960U * sizeof(creal_T));
  c_exp(b_x);
}

/* End of code generation (exp.cpp) */
