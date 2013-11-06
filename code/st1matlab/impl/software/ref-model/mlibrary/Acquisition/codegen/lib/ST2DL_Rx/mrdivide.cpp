/*
 * mrdivide.cpp
 *
 * Code generation for function 'mrdivide'
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
#include "mrdivide.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 *
 */
void b_mrdivide(const real_T A[3960], real_T B, real_T y[3960])
{
  int32_T i10;
  for (i10 = 0; i10 < 3960; i10++) {
    y[i10] = A[i10] / B;
  }
}

/*
 *
 */
creal_T c_mrdivide(const creal_T A, const creal_T B)
{
  creal_T y;
  real_T brm;
  real_T bim;
  real_T d;
  if (B.im == 0.0) {
    if (A.im == 0.0) {
      y.re = A.re / B.re;
      y.im = 0.0;
    } else if (A.re == 0.0) {
      y.re = 0.0;
      y.im = A.im / B.re;
    } else {
      y.re = A.re / B.re;
      y.im = A.im / B.re;
    }
  } else if (B.re == 0.0) {
    if (A.re == 0.0) {
      y.re = A.im / B.im;
      y.im = 0.0;
    } else if (A.im == 0.0) {
      y.re = 0.0;
      y.im = -(A.re / B.im);
    } else {
      y.re = A.im / B.im;
      y.im = -(A.re / B.im);
    }
  } else {
    brm = fabs(B.re);
    bim = fabs(B.im);
    if (brm > bim) {
      bim = B.im / B.re;
      d = B.re + bim * B.im;
      y.re = (A.re + bim * A.im) / d;
      y.im = (A.im - bim * A.re) / d;
    } else if (bim == brm) {
      bim = B.re > 0.0 ? 0.5 : -0.5;
      d = B.im > 0.0 ? 0.5 : -0.5;
      y.re = (A.re * bim + A.im * d) / brm;
      y.im = (A.im * bim - A.re * d) / brm;
    } else {
      bim = B.re / B.im;
      d = B.im + bim * B.re;
      y.re = (bim * A.re + A.im) / d;
      y.im = (bim * A.im - A.re) / d;
    }
  }

  return y;
}

/*
 *
 */
real_T mrdivide(real_T A, real_T B)
{
  return A / B;
}

/* End of code generation (mrdivide.cpp) */
