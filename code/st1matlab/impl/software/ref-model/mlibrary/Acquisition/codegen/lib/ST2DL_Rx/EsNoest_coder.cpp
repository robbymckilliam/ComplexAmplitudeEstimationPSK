/*
 * EsNoest_coder.cpp
 *
 * Code generation for function 'EsNoest_coder'
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
#include "abs.h"
#include "mrdivide.h"
#include "std.h"
#include "mean.h"
#include "exp.h"
#include "mod.h"
#include "round.h"
#include "angle.h"
#include "ST2DL_Rx_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 * function [EsNohat] = EsNoest_coder(MFoutput,M,D,P,p)
 */
real_T EsNoest_coder(const creal_T MFoutput[3990], real_T M, const real_T D[3960],
                     const real_T P[30], const creal_T p[30])
{
  static creal_T symbhat[3990];
  int32_T i9;
  static creal_T b_MFoutput[3960];
  real_T symbsectorhat[3960];
  real_T b_symbsectorhat[3960];
  real_T c_symbsectorhat[3960];
  static const creal_T dc2 = { 0.0, 6.2831853071795862 };

  creal_T y;
  creal_T dcv0[3960];
  real_T symbhat_re;
  real_T symbhat_im;

  /*  Joint pilot and data based Es/No estimator. */
  /* ------------------------------------------------------------- */
  /*  EsNoest.m  */
  /*  */
  /*  [EsNohat] = EsNoest(MFoutput,AcqParams) */
  /*  */
  /*  This function etimates the Es/No from the complex matched filter (MF) */
  /*  output. First, hard decisions for the data symbols are obtained from the */
  /*  MF output. After inserting the known pilot symbols, the symbol stream is */
  /*  used to strip the modulation of the MF output. The Es/No is estimated */
  /*  from the resulting signal. */
  /*   */
  /*  Inputs: */
  /*    MFoutput    - Vector containing complex matched filter output at one  */
  /*                  sample per symbol. */
  /*    AcqParams   - Struct holding aquisition parameters */
  /*        .P              - position (indices) of pilot symbols */
  /*        .p              - pilot symbols */
  /*        .D              - position (indices) of data symbols */
  /*        .M              - M for M-PSK, i.e. size of constellation */
  /*  */
  /*  Outputs: */
  /*    EsNohat     - Estimate of the Es/No. */
  /*  */
  /*  Comments: */
  /*   */
  /* ------------------------------------------------------------- */
  /*  Author: A. Pollok, Institute for Telecommunications Research */
  /*  Project: ASRP */
  /* ------------------------------------------------------------- */
  /*  Copyright 2012  */
  /*  Institute for Telecommunications Research */
  /*  University of South Australia */
  /* ------------------------------------------------------------- */
  /* 'EsNoest_coder:38' symbhat = zeros(1,length(P)+length(D))+1i*zeros(1,length(P)+length(D)); */
  for (i9 = 0; i9 < 3990; i9++) {
    symbhat[i9].re = 0.0;
    symbhat[i9].im = 0.0;
  }

  /*  make hard decisions for PSK data symbols */
  /* 'EsNoest_coder:40' symbsectorhat = mod(round(angle(MFoutput(D))/(2*pi/M)), M); */
  for (i9 = 0; i9 < 3960; i9++) {
    b_MFoutput[i9] = MFoutput[(int32_T)D[i9] - 1];
  }

  angle(b_MFoutput, symbsectorhat);
  memcpy((void *)&b_symbsectorhat[0], (void *)&symbsectorhat[0], 3960U * sizeof
         (real_T));
  b_mrdivide(b_symbsectorhat, mrdivide(6.2831853071795862, M), symbsectorhat);
  b_round(symbsectorhat);
  memcpy((void *)&c_symbsectorhat[0], (void *)&symbsectorhat[0], 3960U * sizeof
         (real_T));
  b_mod(c_symbsectorhat, M, symbsectorhat);

  /* 'EsNoest_coder:41' symbhat(D) = exp(2i*pi/M*symbsectorhat); */
  y = eml_div(dc2, M);
  for (i9 = 0; i9 < 3960; i9++) {
    b_MFoutput[i9].re = symbsectorhat[i9] * y.re;
    b_MFoutput[i9].im = symbsectorhat[i9] * y.im;
  }

  b_exp(b_MFoutput, dcv0);
  for (i9 = 0; i9 < 3960; i9++) {
    symbhat[(int32_T)D[i9] - 1] = dcv0[i9];
  }

  /*  insert pilot symbols */
  /* 'EsNoest_coder:43' symbhat(P) = p; */
  for (i9 = 0; i9 < 30; i9++) {
    symbhat[(int32_T)P[i9] - 1] = p[i9];
  }

  /*  Bill's blind Es/No estimator */
  /* 'EsNoest_coder:46' MFoutDemod = MFoutput .* conj(symbhat); */
  for (i9 = 0; i9 < 3990; i9++) {
    symbhat_re = symbhat[i9].re;
    symbhat_im = -symbhat[i9].im;
    symbhat[i9].re = MFoutput[i9].re * symbhat[i9].re - MFoutput[i9].im *
      -symbhat[i9].im;
    symbhat[i9].im = MFoutput[i9].re * symbhat_im + MFoutput[i9].im * symbhat_re;
  }

  /*  strip modulation off */
  /* 'EsNoest_coder:47' EsNohat = abs(mean(MFoutDemod) / std(MFoutDemod))^2; */
  return rt_powd_snf(b_abs(c_mrdivide(mean(symbhat), b_std(symbhat))), 2.0);
}

/* End of code generation (EsNoest_coder.cpp) */
