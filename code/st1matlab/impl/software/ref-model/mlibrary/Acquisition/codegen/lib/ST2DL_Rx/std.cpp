/*
 * std.cpp
 *
 * Code generation for function 'std'
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
#include "std.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 *
 */
creal_T b_std(const creal_T varargin_1[3990])
{
  creal_T y;
  int32_T ix;
  creal_T xbar;
  int32_T k;
  real_T c_re;
  real_T c_im;
  real_T b_y;
  ix = 0;
  xbar = varargin_1[0];
  for (k = 0; k < 3989; k++) {
    ix++;
    xbar.re += varargin_1[ix].re;
    xbar.im += varargin_1[ix].im;
  }

  xbar = eml_div(xbar, 3990.0);
  ix = 0;
  c_re = varargin_1[0].re - xbar.re;
  c_im = varargin_1[0].im - xbar.im;
  b_y = c_re * c_re + c_im * c_im;
  for (k = 0; k < 3989; k++) {
    ix++;
    c_re = varargin_1[ix].re - xbar.re;
    c_im = varargin_1[ix].im - xbar.im;
    b_y += c_re * c_re + c_im * c_im;
  }

  b_y /= 3989.0;
  y.re = b_y;
  y.im = 0.0;
  if (y.re < 0.0) {
    c_re = 0.0;
    c_im = sqrt(fabs(y.re));
  } else {
    c_re = sqrt(y.re);
    c_im = 0.0;
  }

  y.re = c_re;
  y.im = c_im;
  return y;
}

/* End of code generation (std.cpp) */
