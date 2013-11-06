/*
 * freqest_periodogram_coder.cpp
 *
 * Code generation for function 'freqest_periodogram_coder'
 *
 * C source code generated on: Thu Nov  8 13:52:45 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dcblockST2DL.h"
#include "freqest_periodogram_coder.h"
#include "mf_coder.h"
#include "dcblockST2DL_emxutil.h"
#include "sum.h"
#include "power.h"
#include "ipermute.h"
#include "dcblockST2DL_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void eml_fft(const emxArray_creal_T *x, uint32_T n, emxArray_creal_T *y);
static real_T rt_atan2d_snf(real_T u0, real_T u1);

/* Function Definitions */
static void eml_fft(const emxArray_creal_T *x, uint32_T n, emxArray_creal_T *y)
{
  int32_T iDelta2;
  int32_T nd2;
  uint32_T ju;
  int32_T nRows;
  emxArray_real_T *costab1q;
  int32_T nRowsM1;
  int32_T ixDelta;
  int32_T nRowsD2;
  int32_T nRowsD4;
  int32_T lastChan;
  real_T e;
  int32_T k;
  emxArray_real_T *costab;
  emxArray_real_T *sintab;
  int32_T b_n;
  int32_T n2;
  int32_T ix;
  int32_T chanStart;
  int32_T i;
  uint32_T c_n;
  boolean_T tst;
  real_T temp_re;
  real_T temp_im;
  int32_T iheight;
  int32_T ihi;
  real_T twid_im;
  iDelta2 = y->size[0];
  y->size[0] = (int32_T)n;
  emxEnsureCapacity((emxArray__common *)y, iDelta2, (int32_T)sizeof(creal_T));
  if (n > (uint32_T)x->size[0]) {
    nd2 = y->size[0];
    iDelta2 = y->size[0];
    y->size[0] = nd2;
    emxEnsureCapacity((emxArray__common *)y, iDelta2, (int32_T)sizeof(creal_T));
    nd2--;
    for (iDelta2 = 0; iDelta2 <= nd2; iDelta2++) {
      y->data[iDelta2].re = 0.0;
      y->data[iDelta2].im = 0.0;
    }
  }

  if (x->size[0] == 0) {
  } else {
    ju = n;
    if (ju > 2147483647U) {
      ju = 2147483647U;
    }

    nRows = (int32_T)ju;
    if (x->size[0] > nRows) {
      nd2 = nRows;
    } else {
      nd2 = x->size[0];
    }

    emxInit_real_T(&costab1q, 2);
    nRowsM1 = nd2 - 1;
    nd2 = x->size[0] - nRowsM1;
    ixDelta = 1 >= nd2 ? 1 : nd2;
    nRowsD2 = nRows / 2;
    nRowsD4 = nRowsD2 / 2;
    lastChan = nRows * (x->size[0] / x->size[0] - 1);
    e = 6.2831853071795862 / (real_T)nRows;
    iDelta2 = costab1q->size[0] * costab1q->size[1];
    costab1q->size[0] = 1;
    costab1q->size[1] = nRowsD4 + 1;
    emxEnsureCapacity((emxArray__common *)costab1q, iDelta2, (int32_T)sizeof
                      (real_T));
    costab1q->data[0] = 1.0;
    nd2 = nRowsD4 / 2;
    for (k = 1; k <= nd2; k++) {
      costab1q->data[k] = cos(e * (real_T)k);
    }

    for (k = nd2 + 1; k <= nRowsD4 - 1; k++) {
      costab1q->data[k] = sin(e * (real_T)(nRowsD4 - k));
    }

    emxInit_real_T(&costab, 2);
    emxInit_real_T(&sintab, 2);
    costab1q->data[nRowsD4] = 0.0;
    b_n = costab1q->size[1] - 1;
    n2 = b_n << 1;
    nd2 = n2 + 1;
    iDelta2 = costab->size[0] * costab->size[1];
    costab->size[0] = 1;
    costab->size[1] = nd2;
    emxEnsureCapacity((emxArray__common *)costab, iDelta2, (int32_T)sizeof
                      (real_T));
    iDelta2 = sintab->size[0] * sintab->size[1];
    sintab->size[0] = 1;
    sintab->size[1] = nd2;
    emxEnsureCapacity((emxArray__common *)sintab, iDelta2, (int32_T)sizeof
                      (real_T));
    costab->data[0] = 1.0;
    sintab->data[0] = 0.0;
    for (k = 1; k <= b_n; k++) {
      costab->data[k] = costab1q->data[k];
      sintab->data[k] = -costab1q->data[b_n - k];
    }

    for (k = b_n + 1; k <= n2; k++) {
      costab->data[k] = -costab1q->data[n2 - k];
      sintab->data[k] = -costab1q->data[k - b_n];
    }

    emxFree_real_T(&costab1q);
    ix = 0;
    chanStart = 0;
    while ((nRows > 0) && (chanStart <= lastChan)) {
      ju = 0U;
      nd2 = chanStart;
      for (i = 1; i <= nRowsM1; i++) {
        y->data[nd2] = x->data[ix];
        c_n = n;
        tst = TRUE;
        while (tst) {
          c_n >>= 1U;
          ju ^= c_n;
          tst = ((int32_T)(ju & c_n) == 0);
        }

        nd2 = chanStart + (int32_T)ju;
        ix++;
      }

      y->data[nd2] = x->data[ix];
      ix += ixDelta;
      nd2 = (chanStart + nRows) - 2;
      if (nRows > 1) {
        for (i = chanStart; i <= nd2; i += 2) {
          temp_re = y->data[i + 1].re;
          temp_im = y->data[i + 1].im;
          y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
          y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
          y->data[i].re += temp_re;
          y->data[i].im += temp_im;
        }
      }

      n2 = 2;
      iDelta2 = 4;
      k = nRowsD4;
      iheight = 1 + ((nRowsD4 - 1) << 2);
      while (k > 0) {
        i = chanStart;
        ihi = chanStart + iheight;
        while (i < ihi) {
          nd2 = i + n2;
          temp_re = y->data[nd2].re;
          temp_im = y->data[nd2].im;
          y->data[i + n2].re = y->data[i].re - y->data[nd2].re;
          y->data[i + n2].im = y->data[i].im - y->data[nd2].im;
          y->data[i].re += temp_re;
          y->data[i].im += temp_im;
          i += iDelta2;
        }

        nd2 = chanStart + 1;
        for (b_n = k; b_n < nRowsD2; b_n += k) {
          e = costab->data[b_n];
          twid_im = sintab->data[b_n];
          i = nd2;
          ihi = nd2 + iheight;
          while (i < ihi) {
            temp_re = e * y->data[i + n2].re - twid_im * y->data[i + n2].im;
            temp_im = e * y->data[i + n2].im + twid_im * y->data[i + n2].re;
            y->data[i + n2].re = y->data[i].re - temp_re;
            y->data[i + n2].im = y->data[i].im - temp_im;
            y->data[i].re += temp_re;
            y->data[i].im += temp_im;
            i += iDelta2;
          }

          nd2++;
        }

        k /= 2;
        n2 = iDelta2;
        iDelta2 <<= 1;
        iheight -= n2;
      }

      chanStart += nRows;
    }

    emxFree_real_T(&sintab);
    emxFree_real_T(&costab);
  }
}

static real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    y = atan2(u0 > 0.0 ? 1.0 : -1.0, u1 > 0.0 ? 1.0 : -1.0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

void freqest_periodogram_coder(const creal_T RxSamp[16000], const real_T tSamp
  [16000], real_T T, real_T OS, real_T Lsin, real_T taumin, real_T taumax,
  real_T *fhat, real_T *phihat)
{
  real_T Ts;
  real_T Yd_re;
  real_T idxEnd;
  int32_T i0;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  emxArray_int32_T *r0;
  int32_T i4;
  int32_T absn;
  emxArray_creal_T *b_RxSamp;
  int32_T iv0[2];
  int32_T ix;
  int32_T itmp;
  emxArray_int32_T r1;
  int32_T n;
  emxArray_creal_T *c_RxSamp;
  real_T Nfft;
  emxArray_creal_T *I;
  emxArray_creal_T *b_y1;
  uint32_T u0;
  uint32_T uv0[2];
  emxArray_real_T *y;
  real_T a;
  real_T b;
  emxArray_real_T *b_y;
  boolean_T exitg1;
  real_T Y_re;
  real_T Y_im;
  emxArray_creal_T *r2;
  emxArray_creal_T *c_y;
  emxArray_real_T *b_tSamp;
  emxArray_creal_T *c_tSamp;
  emxArray_int32_T *r3;
  emxArray_int32_T *r4;
  emxArray_int32_T *r5;
  emxArray_int32_T *r6;
  emxArray_int32_T *r7;
  emxArray_int32_T *r8;
  emxArray_int32_T *r9;
  emxArray_creal_T *d_RxSamp;
  emxArray_creal_T *e_RxSamp;
  emxArray_creal_T *f_RxSamp;
  creal_T b_b;
  int32_T iv1[2];
  int32_T iv2[2];
  real_T Id;
  int32_T iv3[2];

  /*  Frequency estimation from periodogram. */
  /* ------------------------------------------------------------- */
  /*  freqest_periodogram.m  */
  /*  */
  /*  [fhat] = freqest_periodogram(RxSamp,tSamp,AcqParams) */
  /*  */
  /*  This function etimates the frequency of a noisy complex exponential by */
  /*  maximising the periodogram. First, a coarse estimate is obtained by */
  /*  computing the periodogram at the discrete frequencies of an FFT. This */
  /*  coarse estimate then serves as the initialisation for the second */
  /*  estimator stage, which climbs the objective function using Newton's */
  /*  method. As proposed by Quinn et al. [1], Newton's method is applied to */
  /*  the logarithm of the periodogram rather than the periodogram itself. This */
  /*  approach eliminates the need for zero padding to increase the resolution */
  /*  of the FFT-based initial estimate. */
  /*  */
  /*    [1] B. Quinn et al., "Maximizing the Periodogram", in Proc. IEEE Global */
  /*        Global Telecommunications Conference, Nov. 2008  */
  /*   */
  /*  Inputs: */
  /*    RxSamp      - Vector containing samples of received signal */
  /*    tSamp       - Vector of sample times in seconds */
  /*    AcqParams   - Struct holding aquisition parameters */
  /*        .T          - symbol period */
  /*        .OS         - oversampling factor */
  /*        .Lsin       - duration of the sinusoid in symbol periods */
  /*        .taumin         - minimum time offset */
  /*        .taumax         - maximum time offset */
  /*  */
  /*  Outputs: */
  /*   fhat         - Frequency estimate [Hz] */
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
  /* % system parameters */
  /* T = AcqParams.T; % symbol period */
  /* OS = AcqParams.OS; % oversampling factor */
  Ts = T / OS;

  /*  sample period */
  /* Lsin = AcqParams.Lsin; % duration of the sinusoid in symbol periods */
  /*  channel parameters */
  /* taumin = AcqParams.taumin; */
  /* taumax = AcqParams.taumax; */
  /*  discard samples before minimum time delay */
  Yd_re = floor(taumin / Ts);
  idxEnd = (((Yd_re + 1.0) - 1.0) + ceil((taumax - taumin) / Ts)) + Lsin * OS;
  if (Yd_re + 1.0 > idxEnd) {
    i0 = 1;
    i1 = 0;
  } else {
    i0 = (int32_T)(Yd_re + 1.0);
    i1 = (int32_T)idxEnd;
  }

  if (Yd_re + 1.0 > idxEnd) {
    i2 = 1;
    i3 = 0;
  } else {
    i2 = (int32_T)(Yd_re + 1.0);
    i3 = (int32_T)idxEnd;
  }

  emxInit_int32_T(&r0, 1);

  /*  number of samples */
  /* % coarse frequency estimation via periodogram */
  /*  zero-padding factor for FFT */
  i4 = r0->size[0];
  r0->size[0] = (i1 - i0) + 1;
  emxEnsureCapacity((emxArray__common *)r0, i4, (int32_T)sizeof(int32_T));
  absn = i1 - i0;
  for (i4 = 0; i4 <= absn; i4++) {
    r0->data[i4] = i0 + i4;
  }

  emxInit_creal_T(&b_RxSamp, 2);
  iv0[0] = 1;
  iv0[1] = r0->size[0];
  i4 = b_RxSamp->size[0] * b_RxSamp->size[1];
  b_RxSamp->size[0] = iv0[0];
  b_RxSamp->size[1] = iv0[1];
  emxEnsureCapacity((emxArray__common *)b_RxSamp, i4, (int32_T)sizeof(creal_T));
  absn = iv0[1] - 1;
  for (i4 = 0; i4 <= absn; i4++) {
    ix = iv0[0] - 1;
    for (itmp = 0; itmp <= ix; itmp++) {
      r1 = *r0;
      r1.size = (int32_T *)&iv0;
      r1.numDimensions = 1;
      b_RxSamp->data[itmp + b_RxSamp->size[0] * i4] = RxSamp[r1.data[itmp +
        r1.size[0] * i4] - 1];
    }
  }

  emxFree_int32_T(&r0);
  absn = b_RxSamp->size[1];
  Yd_re = frexp((real_T)absn, &n);
  idxEnd = (real_T)n;
  emxFree_creal_T(&b_RxSamp);
  if (Yd_re == 0.5) {
    idxEnd = (real_T)n - 1.0;
  }

  b_emxInit_creal_T(&c_RxSamp, 1);
  Nfft = rt_powd_snf(2.0, idxEnd);

  /*  compute periodogram */
  i4 = c_RxSamp->size[0];
  c_RxSamp->size[0] = (i1 - i0) + 1;
  emxEnsureCapacity((emxArray__common *)c_RxSamp, i4, (int32_T)sizeof(creal_T));
  absn = i1 - i0;
  for (i4 = 0; i4 <= absn; i4++) {
    c_RxSamp->data[i4] = RxSamp[(i0 + i4) - 1];
  }

  emxInit_creal_T(&I, 2);
  b_emxInit_creal_T(&b_y1, 1);
  Yd_re = Nfft;
  Yd_re = floor(Yd_re + 0.5);
  if (Yd_re < 4.294967296E+9) {
    u0 = (uint32_T)Yd_re;
  } else {
    u0 = MAX_uint32_T;
  }

  eml_fft(c_RxSamp, u0, b_y1);
  ipermute(b_y1, I);
  emxFree_creal_T(&c_RxSamp);
  emxFree_creal_T(&b_y1);
  for (i4 = 0; i4 < 2; i4++) {
    uv0[i4] = (uint32_T)I->size[i4];
  }

  emxInit_real_T(&y, 2);
  i4 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int32_T)uv0[1];
  emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
  for (absn = 0; absn <= I->size[1] - 1; absn++) {
    a = fabs(I->data[(int32_T)(1.0 + (real_T)absn) - 1].re);
    b = fabs(I->data[(int32_T)(1.0 + (real_T)absn) - 1].im);
    if (a < b) {
      a /= b;
      b *= sqrt(a * a + 1.0);
    } else if (a > b) {
      b /= a;
      b = sqrt(b * b + 1.0) * a;
    } else if (rtIsNaN(b)) {
    } else {
      b = a * 1.4142135623730951;
    }

    y->data[(int32_T)(1.0 + (real_T)absn) - 1] = b;
  }

  emxInit_real_T(&b_y, 2);
  i4 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
  absn = y->size[0] * y->size[1] - 1;
  for (i4 = 0; i4 <= absn; i4++) {
    b_y->data[i4] = y->data[i4];
  }

  power(b_y, y);
  absn = 1;
  n = y->size[1];
  Yd_re = y->data[0];
  itmp = 0;
  emxFree_real_T(&b_y);
  if (n > 1) {
    if (rtIsNaN(y->data[0])) {
      ix = 2;
      exitg1 = 0U;
      while ((exitg1 == 0U) && (ix <= n)) {
        absn = ix;
        if (!rtIsNaN(y->data[ix - 1])) {
          Yd_re = y->data[ix - 1];
          exitg1 = 1U;
        } else {
          ix++;
        }
      }
    }

    if (absn < n) {
      while (absn + 1 <= n) {
        if (y->data[absn] > Yd_re) {
          Yd_re = y->data[absn];
          itmp = absn;
        }

        absn++;
      }
    }
  }

  /* % refine estimate via Newton's method */
  /*  alpha=0: max log(periodogram), alpha=1: max standard periodogram */
  idxEnd = rtInf;
  n = 0;
  *fhat = ((real_T)(itmp + 1) - 1.0) / Nfft / Ts;
  Y_re = 0.0;
  Y_im = 0.0;
  emxInit_creal_T(&r2, 2);
  emxInit_creal_T(&c_y, 2);
  emxInit_real_T(&b_tSamp, 2);
  emxInit_creal_T(&c_tSamp, 2);
  emxInit_int32_T(&r3, 1);
  emxInit_int32_T(&r4, 1);
  emxInit_int32_T(&r5, 1);
  emxInit_int32_T(&r6, 1);
  emxInit_int32_T(&r7, 1);
  emxInit_int32_T(&r8, 1);
  emxInit_int32_T(&r9, 1);
  emxInit_creal_T(&d_RxSamp, 2);
  emxInit_creal_T(&e_RxSamp, 2);
  emxInit_creal_T(&f_RxSamp, 2);
  while ((fabs(idxEnd) > 1.0E-10) && (n < 15)) {
    /*  evaluate periodogram I and its first two derivatives Idd and Idd */
    Y_re = *fhat * -0.0;
    Y_im = *fhat * -6.2831853071795862;
    i4 = r3->size[0];
    r3->size[0] = (i3 - i2) + 1;
    emxEnsureCapacity((emxArray__common *)r3, i4, (int32_T)sizeof(int32_T));
    absn = i3 - i2;
    for (i4 = 0; i4 <= absn; i4++) {
      r3->data[i4] = i2 + i4;
    }

    i4 = I->size[0] * I->size[1];
    I->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)I, i4, (int32_T)sizeof(creal_T));
    absn = r3->size[0];
    i4 = I->size[0] * I->size[1];
    I->size[1] = absn;
    emxEnsureCapacity((emxArray__common *)I, i4, (int32_T)sizeof(creal_T));
    absn = r3->size[0] - 1;
    for (i4 = 0; i4 <= absn; i4++) {
      I->data[i4].re = tSamp[r3->data[i4] - 1] * Y_re;
      I->data[i4].im = tSamp[r3->data[i4] - 1] * Y_im;
    }

    i4 = r2->size[0] * r2->size[1];
    r2->size[0] = 1;
    r2->size[1] = I->size[1];
    emxEnsureCapacity((emxArray__common *)r2, i4, (int32_T)sizeof(creal_T));
    absn = I->size[0] * I->size[1] - 1;
    for (i4 = 0; i4 <= absn; i4++) {
      r2->data[i4] = I->data[i4];
    }

    for (absn = 0; absn <= I->size[1] - 1; absn++) {
      Yd_re = exp(r2->data[absn].re / 2.0);
      a = r2->data[absn].im;
      idxEnd = r2->data[absn].im;
      r2->data[absn].re = Yd_re * (Yd_re * cos(a));
      r2->data[absn].im = Yd_re * (Yd_re * sin(idxEnd));
    }

    i4 = r4->size[0];
    r4->size[0] = (i1 - i0) + 1;
    emxEnsureCapacity((emxArray__common *)r4, i4, (int32_T)sizeof(int32_T));
    absn = i1 - i0;
    for (i4 = 0; i4 <= absn; i4++) {
      r4->data[i4] = i0 + i4;
    }

    i4 = I->size[0] * I->size[1];
    I->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)I, i4, (int32_T)sizeof(creal_T));
    absn = r4->size[0];
    i4 = I->size[0] * I->size[1];
    I->size[1] = absn;
    emxEnsureCapacity((emxArray__common *)I, i4, (int32_T)sizeof(creal_T));
    absn = r4->size[0] - 1;
    for (i4 = 0; i4 <= absn; i4++) {
      idxEnd = RxSamp[r4->data[i4] - 1].re;
      Yd_re = RxSamp[r4->data[i4] - 1].im;
      Nfft = r2->data[i4].re;
      a = r2->data[i4].im;
      I->data[i4].re = idxEnd * Nfft - Yd_re * a;
      I->data[i4].im = idxEnd * a + Yd_re * Nfft;
    }

    b_b = sum(I);
    i4 = r5->size[0];
    r5->size[0] = (i1 - i0) + 1;
    emxEnsureCapacity((emxArray__common *)r5, i4, (int32_T)sizeof(int32_T));
    absn = i1 - i0;
    for (i4 = 0; i4 <= absn; i4++) {
      r5->data[i4] = i0 + i4;
    }

    iv1[0] = 1;
    iv1[1] = r5->size[0];
    i4 = d_RxSamp->size[0] * d_RxSamp->size[1];
    d_RxSamp->size[0] = iv1[0];
    d_RxSamp->size[1] = iv1[1];
    emxEnsureCapacity((emxArray__common *)d_RxSamp, i4, (int32_T)sizeof(creal_T));
    absn = iv1[1] - 1;
    for (i4 = 0; i4 <= absn; i4++) {
      ix = iv1[0] - 1;
      for (itmp = 0; itmp <= ix; itmp++) {
        r1 = *r5;
        r1.size = (int32_T *)&iv1;
        r1.numDimensions = 1;
        d_RxSamp->data[itmp + d_RxSamp->size[0] * i4] = RxSamp[r1.data[itmp +
          r1.size[0] * i4] - 1];
      }
    }

    Y_re = 1.0 / (real_T)d_RxSamp->size[1] * b_b.re;
    Y_im = 1.0 / (real_T)d_RxSamp->size[1] * b_b.im;
    i4 = r6->size[0];
    r6->size[0] = (i3 - i2) + 1;
    emxEnsureCapacity((emxArray__common *)r6, i4, (int32_T)sizeof(int32_T));
    absn = i3 - i2;
    for (i4 = 0; i4 <= absn; i4++) {
      r6->data[i4] = i2 + i4;
    }

    i4 = c_tSamp->size[0] * c_tSamp->size[1];
    c_tSamp->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)c_tSamp, i4, (int32_T)sizeof(creal_T));
    absn = r6->size[0];
    i4 = c_tSamp->size[0] * c_tSamp->size[1];
    c_tSamp->size[1] = absn;
    emxEnsureCapacity((emxArray__common *)c_tSamp, i4, (int32_T)sizeof(creal_T));
    absn = r6->size[0] - 1;
    for (i4 = 0; i4 <= absn; i4++) {
      c_tSamp->data[i4].re = tSamp[r6->data[i4] - 1] * I->data[i4].re;
      c_tSamp->data[i4].im = tSamp[r6->data[i4] - 1] * I->data[i4].im;
    }

    b_b = sum(c_tSamp);
    i4 = r7->size[0];
    r7->size[0] = (i1 - i0) + 1;
    emxEnsureCapacity((emxArray__common *)r7, i4, (int32_T)sizeof(int32_T));
    absn = i1 - i0;
    for (i4 = 0; i4 <= absn; i4++) {
      r7->data[i4] = i0 + i4;
    }

    iv2[0] = 1;
    iv2[1] = r7->size[0];
    i4 = e_RxSamp->size[0] * e_RxSamp->size[1];
    e_RxSamp->size[0] = iv2[0];
    e_RxSamp->size[1] = iv2[1];
    emxEnsureCapacity((emxArray__common *)e_RxSamp, i4, (int32_T)sizeof(creal_T));
    absn = iv2[1] - 1;
    for (i4 = 0; i4 <= absn; i4++) {
      ix = iv2[0] - 1;
      for (itmp = 0; itmp <= ix; itmp++) {
        r1 = *r7;
        r1.size = (int32_T *)&iv2;
        r1.numDimensions = 1;
        e_RxSamp->data[itmp + e_RxSamp->size[0] * i4] = RxSamp[r1.data[itmp +
          r1.size[0] * i4] - 1];
      }
    }

    absn = e_RxSamp->size[1];
    a = -6.2831853071795862 / (real_T)absn;
    Yd_re = 0.0 * b_b.re - a * b_b.im;
    idxEnd = 0.0 * b_b.im + a * b_b.re;
    i4 = r8->size[0];
    r8->size[0] = (i3 - i2) + 1;
    emxEnsureCapacity((emxArray__common *)r8, i4, (int32_T)sizeof(int32_T));
    absn = i3 - i2;
    for (i4 = 0; i4 <= absn; i4++) {
      r8->data[i4] = i2 + i4;
    }

    i4 = b_tSamp->size[0] * b_tSamp->size[1];
    b_tSamp->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)b_tSamp, i4, (int32_T)sizeof(real_T));
    absn = r8->size[0];
    i4 = b_tSamp->size[0] * b_tSamp->size[1];
    b_tSamp->size[1] = absn;
    emxEnsureCapacity((emxArray__common *)b_tSamp, i4, (int32_T)sizeof(real_T));
    absn = r8->size[0] - 1;
    for (i4 = 0; i4 <= absn; i4++) {
      b_tSamp->data[i4] = tSamp[r8->data[i4] - 1];
    }

    power(b_tSamp, y);
    i4 = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    c_y->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(creal_T));
    absn = y->size[0] * y->size[1] - 1;
    for (i4 = 0; i4 <= absn; i4++) {
      c_y->data[i4].re = y->data[i4] * I->data[i4].re;
      c_y->data[i4].im = y->data[i4] * I->data[i4].im;
    }

    b_b = sum(c_y);
    a = fabs(Y_re);
    b = fabs(Y_im);
    if (a < b) {
      a /= b;
      b *= sqrt(a * a + 1.0);
    } else if (a > b) {
      b /= a;
      b = sqrt(b * b + 1.0) * a;
    } else if (rtIsNaN(b)) {
    } else {
      b = a * 1.4142135623730951;
    }

    Id = 2.0 * (Yd_re * Y_re - idxEnd * -Y_im);
    a = fabs(Yd_re);
    idxEnd = fabs(idxEnd);
    if (a < idxEnd) {
      a /= idxEnd;
      idxEnd *= sqrt(a * a + 1.0);
    } else if (a > idxEnd) {
      idxEnd /= a;
      idxEnd = sqrt(idxEnd * idxEnd + 1.0) * a;
    } else if (rtIsNaN(idxEnd)) {
    } else {
      idxEnd = a * 1.4142135623730951;
    }

    /*  Newton iterate */
    i4 = r9->size[0];
    r9->size[0] = (i1 - i0) + 1;
    emxEnsureCapacity((emxArray__common *)r9, i4, (int32_T)sizeof(int32_T));
    absn = i1 - i0;
    for (i4 = 0; i4 <= absn; i4++) {
      r9->data[i4] = i0 + i4;
    }

    iv3[0] = 1;
    iv3[1] = r9->size[0];
    i4 = f_RxSamp->size[0] * f_RxSamp->size[1];
    f_RxSamp->size[0] = iv3[0];
    f_RxSamp->size[1] = iv3[1];
    emxEnsureCapacity((emxArray__common *)f_RxSamp, i4, (int32_T)sizeof(creal_T));
    absn = iv3[1] - 1;
    for (i4 = 0; i4 <= absn; i4++) {
      ix = iv3[0] - 1;
      for (itmp = 0; itmp <= ix; itmp++) {
        r1 = *r9;
        r1.size = (int32_T *)&iv3;
        r1.numDimensions = 1;
        f_RxSamp->data[itmp + f_RxSamp->size[0] * i4] = RxSamp[r1.data[itmp +
          r1.size[0] * i4] - 1];
      }
    }

    Nfft = -39.478417604357432 / (real_T)f_RxSamp->size[1] * b_b.re;
    a = -39.478417604357432 / (real_T)f_RxSamp->size[1] * b_b.im;
    idxEnd = Id / ((2.0 * (Nfft * Y_re - a * -Y_im) + 2.0 * rt_powd_snf(idxEnd,
      2.0)) + -rt_powd_snf(Id, 2.0) / rt_powd_snf(b, 2.0));
    *fhat -= idxEnd;
    n++;
  }

  emxFree_creal_T(&f_RxSamp);
  emxFree_creal_T(&e_RxSamp);
  emxFree_creal_T(&d_RxSamp);
  emxFree_int32_T(&r9);
  emxFree_int32_T(&r8);
  emxFree_int32_T(&r7);
  emxFree_int32_T(&r6);
  emxFree_int32_T(&r5);
  emxFree_int32_T(&r4);
  emxFree_int32_T(&r3);
  emxFree_creal_T(&c_tSamp);
  emxFree_real_T(&b_tSamp);
  emxFree_creal_T(&c_y);
  emxFree_real_T(&y);
  emxFree_creal_T(&r2);
  emxFree_creal_T(&I);

  /* % handle negative frequency offsets */
  if (*fhat > 1.0 / (2.0 * Ts)) {
    /*  negative frequency offset */
    *fhat -= 1.0 / Ts;
  }

  /* % estimate complex amplitude (least squares estimator) */
  /*  amplsin = acqparams.amplsin; */
  /*  cgainhat = Y/amplsin; */
  *phihat = rt_atan2d_snf(Y_im, Y_re);
}

/* End of code generation (freqest_periodogram_coder.cpp) */
