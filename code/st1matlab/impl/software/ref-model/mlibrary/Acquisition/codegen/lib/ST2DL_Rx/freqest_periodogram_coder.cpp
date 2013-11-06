/*
 * freqest_periodogram_coder.cpp
 *
 * Code generation for function 'freqest_periodogram_coder'
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
#include "ST2DL_Rx_emxutil.h"
#include "sum.h"
#include "power.h"
#include "ipermute.h"
#include "ST2DL_Rx_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void eml_fft(const emxArray_creal_T *x, uint32_T n, emxArray_creal_T *y);

/* Function Definitions */

/*
 *
 */
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

/*
 *
 */
creal_T eml_div(const creal_T x, real_T y)
{
  creal_T z;
  if (x.im == 0.0) {
    z.re = x.re / y;
    z.im = 0.0;
  } else if (x.re == 0.0) {
    z.re = 0.0;
    z.im = x.im / y;
  } else {
    z.re = x.re / y;
    z.im = x.im / y;
  }

  return z;
}

/*
 * function [fhat,phihat] = freqest_periodogram_coder(RxSamp,tSamp,T,OS,Lsin,taumin,taumax)
 */
void freqest_periodogram_coder(const creal_T RxSamp[16000], const real_T tSamp
  [16000], real_T T, real_T OS, real_T Lsin, real_T taumin, real_T taumax,
  real_T *fhat, real_T *phihat)
{
  real_T Ts;
  real_T Yd_im;
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
  static const creal_T dc0 = { -0.0, -6.2831853071795862 };

  creal_T dc1;
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
  /* 'freqest_periodogram_coder:49' Ts = T/OS; */
  Ts = T / OS;

  /*  sample period */
  /* Lsin = AcqParams.Lsin; % duration of the sinusoid in symbol periods */
  /*  channel parameters */
  /* taumin = AcqParams.taumin; */
  /* taumax = AcqParams.taumax; */
  /*  discard samples before minimum time delay */
  /* 'freqest_periodogram_coder:57' idxStart = floor(taumin/Ts)+1; */
  Yd_im = floor(taumin / Ts);

  /* 'freqest_periodogram_coder:58' idxEnd   = idxStart-1 + ceil((taumax-taumin)/Ts) + Lsin*OS; */
  idxEnd = (((Yd_im + 1.0) - 1.0) + ceil((taumax - taumin) / Ts)) + Lsin * OS;

  /* 'freqest_periodogram_coder:59' RxSamp = RxSamp(idxStart:idxEnd); */
  if (Yd_im + 1.0 > idxEnd) {
    i0 = 1;
    i1 = 0;
  } else {
    i0 = (int32_T)(Yd_im + 1.0);
    i1 = (int32_T)idxEnd;
  }

  /* 'freqest_periodogram_coder:60' tSamp  =  tSamp(idxStart:idxEnd); */
  if (Yd_im + 1.0 > idxEnd) {
    i2 = 1;
    i3 = 0;
  } else {
    i2 = (int32_T)(Yd_im + 1.0);
    i3 = (int32_T)idxEnd;
  }

  emxInit_int32_T(&r0, 1);

  /* 'freqest_periodogram_coder:62' N = length(RxSamp); */
  /*  number of samples */
  /* % coarse frequency estimation via periodogram */
  /* 'freqest_periodogram_coder:65' pad = 1; */
  /*  zero-padding factor for FFT */
  /* 'freqest_periodogram_coder:66' Nfft = 2^nextpow2(pad*N); */
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
  Yd_im = frexp((real_T)absn, &n);
  idxEnd = (real_T)n;
  emxFree_creal_T(&b_RxSamp);
  if (Yd_im == 0.5) {
    idxEnd = (real_T)n - 1.0;
  }

  b_emxInit_creal_T(&c_RxSamp, 1);
  Nfft = rt_powd_snf(2.0, idxEnd);

  /*  compute periodogram */
  /* 'freqest_periodogram_coder:68' I = fft(RxSamp,Nfft); */
  i4 = c_RxSamp->size[0];
  c_RxSamp->size[0] = (i1 - i0) + 1;
  emxEnsureCapacity((emxArray__common *)c_RxSamp, i4, (int32_T)sizeof(creal_T));
  absn = i1 - i0;
  for (i4 = 0; i4 <= absn; i4++) {
    c_RxSamp->data[i4] = RxSamp[(i0 + i4) - 1];
  }

  emxInit_creal_T(&I, 2);
  b_emxInit_creal_T(&b_y1, 1);
  Yd_im = Nfft;
  Yd_im = floor(Yd_im + 0.5);
  if (Yd_im < 4.294967296E+9) {
    u0 = (uint32_T)Yd_im;
  } else {
    u0 = MAX_uint32_T;
  }

  eml_fft(c_RxSamp, u0, b_y1);
  ipermute(b_y1, I);

  /* 'freqest_periodogram_coder:69' [~,idxmax] = max(abs(I).^2); */
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
  Yd_im = y->data[0];
  itmp = 0;
  emxFree_real_T(&b_y);
  if (n > 1) {
    if (rtIsNaN(y->data[0])) {
      ix = 2;
      exitg1 = 0U;
      while ((exitg1 == 0U) && (ix <= n)) {
        absn = ix;
        if (!rtIsNaN(y->data[ix - 1])) {
          Yd_im = y->data[ix - 1];
          exitg1 = 1U;
        } else {
          ix++;
        }
      }
    }

    if (absn < n) {
      while (absn + 1 <= n) {
        if (y->data[absn] > Yd_im) {
          Yd_im = y->data[absn];
          itmp = absn;
        }

        absn++;
      }
    }
  }

  /* 'freqest_periodogram_coder:70' fhatcoarse = (idxmax-1)/Nfft/Ts; */
  /* % refine estimate via Newton's method */
  /* 'freqest_periodogram_coder:73' alpha = 0; */
  /*  alpha=0: max log(periodogram), alpha=1: max standard periodogram */
  /* 'freqest_periodogram_coder:74' gamma = Inf; */
  idxEnd = rtInf;

  /* 'freqest_periodogram_coder:75' EPSILON = 1e-10; */
  /* 'freqest_periodogram_coder:76' MAX_ITER = 15; */
  /* 'freqest_periodogram_coder:77' numIter = 0; */
  n = 0;

  /* 'freqest_periodogram_coder:78' fhat = fhatcoarse; */
  *fhat = ((real_T)(itmp + 1) - 1.0) / Nfft / Ts;

  /* 'freqest_periodogram_coder:79' Y=0; */
  Y_re = 0.0;
  Y_im = 0.0;

  /* 'freqest_periodogram_coder:80' while ((abs(gamma) > EPSILON) && (numIter < MAX_ITER)) */
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
    /* 'freqest_periodogram_coder:82' X = RxSamp.*exp(-2i*pi*fhat*tSamp); */
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
      Yd_im = exp(r2->data[absn].re / 2.0);
      a = r2->data[absn].im;
      idxEnd = r2->data[absn].im;
      r2->data[absn].re = Yd_im * (Yd_im * cos(a));
      r2->data[absn].im = Yd_im * (Yd_im * sin(idxEnd));
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
      Yd_im = RxSamp[r4->data[i4] - 1].im;
      Nfft = r2->data[i4].re;
      a = r2->data[i4].im;
      I->data[i4].re = idxEnd * Nfft - Yd_im * a;
      I->data[i4].im = idxEnd * a + Yd_im * Nfft;
    }

    /* 'freqest_periodogram_coder:83' Y = 1/N*sum(X); */
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

    /* 'freqest_periodogram_coder:84' Yd =  -2i*pi/N  * sum(tSamp   .*X); */
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

    dc1 = eml_div(dc0, (real_T)e_RxSamp->size[1]);
    idxEnd = dc1.re * b_b.re - dc1.im * b_b.im;
    Yd_im = dc1.re * b_b.im + dc1.im * b_b.re;

    /* 'freqest_periodogram_coder:85' Ydd = -4*pi^2/N * sum(tSamp.^2.*X); */
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

    /* 'freqest_periodogram_coder:86' I = abs(Y)^2; */
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

    /* 'freqest_periodogram_coder:87' Id =  2*real(Yd*conj(Y)); */
    Id = 2.0 * (idxEnd * Y_re - Yd_im * -Y_im);

    /* 'freqest_periodogram_coder:88' Idd = 2*real(Ydd*conj(Y)) + 2*abs(Yd)^2; */
    a = fabs(idxEnd);
    idxEnd = fabs(Yd_im);
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
    /* 'freqest_periodogram_coder:90' gamma = Id/(Idd + (alpha-1)*Id^2/I); */
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

    /* 'freqest_periodogram_coder:91' fhat = fhat - gamma; */
    *fhat -= idxEnd;

    /* 'freqest_periodogram_coder:92' numIter = numIter + 1; */
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
  /* 'freqest_periodogram_coder:96' if (fhat > 1/(2*Ts)) */
  if (*fhat > 1.0 / (2.0 * Ts)) {
    /*  negative frequency offset */
    /* 'freqest_periodogram_coder:98' fhat = fhat - 1/Ts; */
    *fhat -= 1.0 / Ts;
  }

  /* % estimate complex amplitude (least squares estimator) */
  /*  amplsin = acqparams.amplsin; */
  /*  cgainhat = Y/amplsin; */
  /* 'freqest_periodogram_coder:104' phihat = angle(Y); */
  *phihat = rt_atan2d_snf(Y_im, Y_re);
}

/* End of code generation (freqest_periodogram_coder.cpp) */
