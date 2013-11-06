/*
 * mf_coder.cpp
 *
 * Code generation for function 'mf_coder'
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
#include "ipermute.h"
#include "abs.h"
#include "rdivide.h"
#include "power.h"
#include "mrdivide.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_eml_fft(const emxArray_real_T *x, uint32_T n, emxArray_creal_T *y);
static void c_eml_fft(const creal_T x[16000], uint32_T n, emxArray_creal_T *y);
static void d_eml_fft(const emxArray_creal_T *x, uint32_T n, emxArray_creal_T *y);
static void eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y);

/* Function Definitions */
static void b_eml_fft(const emxArray_real_T *x, uint32_T n, emxArray_creal_T *y)
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
        y->data[nd2].re = x->data[ix];
        y->data[nd2].im = 0.0;
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

      y->data[nd2].re = x->data[ix];
      y->data[nd2].im = 0.0;
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

static void c_eml_fft(const creal_T x[16000], uint32_T n, emxArray_creal_T *y)
{
  uint32_T sz[2];
  int32_T iDelta2;
  int32_T nd2;
  uint32_T ju;
  int32_T nRows;
  int32_T minval;
  emxArray_real_T *costab1q;
  int32_T ixDelta;
  int32_T nRowsD2;
  int32_T nRowsD4;
  real_T e;
  int32_T k;
  emxArray_real_T *costab;
  emxArray_real_T *sintab;
  int32_T b_n;
  int32_T n2;
  int32_T ix;
  int32_T i;
  uint32_T c_n;
  boolean_T tst;
  real_T temp_re;
  real_T temp_im;
  int32_T iheight;
  int32_T ihi;
  real_T twid_im;
  for (iDelta2 = 0; iDelta2 < 2; iDelta2++) {
    sz[iDelta2] = (uint32_T)(16000 + -15999 * iDelta2);
  }

  sz[0] = n;
  iDelta2 = y->size[0];
  y->size[0] = (int32_T)sz[0];
  emxEnsureCapacity((emxArray__common *)y, iDelta2, (int32_T)sizeof(creal_T));
  if (n > 16000U) {
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

  ju = n;
  if (ju > 2147483647U) {
    ju = 2147483647U;
  }

  nRows = (int32_T)ju;
  if (16000 > nRows) {
    minval = nRows;
  } else {
    minval = 16000;
  }

  emxInit_real_T(&costab1q, 2);
  nd2 = 16001 - minval;
  ixDelta = 1 >= nd2 ? 1 : nd2;
  nRowsD2 = nRows / 2;
  nRowsD4 = nRowsD2 / 2;
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
  emxEnsureCapacity((emxArray__common *)costab, iDelta2, (int32_T)sizeof(real_T));
  iDelta2 = sintab->size[0] * sintab->size[1];
  sintab->size[0] = 1;
  sintab->size[1] = nd2;
  emxEnsureCapacity((emxArray__common *)sintab, iDelta2, (int32_T)sizeof(real_T));
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
  nd2 = 0;
  while ((nRows > 0) && (nd2 <= 0)) {
    ju = 0U;
    nd2 = 0;
    for (i = 1; i <= minval - 1; i++) {
      y->data[nd2] = x[ix];
      c_n = n;
      tst = TRUE;
      while (tst) {
        c_n >>= 1U;
        ju ^= c_n;
        tst = ((int32_T)(ju & c_n) == 0);
      }

      nd2 = (int32_T)ju;
      ix++;
    }

    y->data[nd2] = x[ix];
    ix += ixDelta;
    nd2 = nRows;
    if (nRows > 1) {
      for (i = 0; i <= nd2 - 2; i += 2) {
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
      i = 0;
      ihi = iheight;
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

      nd2 = 1;
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

    nd2 = nRows;
  }

  emxFree_real_T(&sintab);
  emxFree_real_T(&costab);
}

static void d_eml_fft(const emxArray_creal_T *x, uint32_T n, emxArray_creal_T *y)
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
      sintab->data[k] = costab1q->data[b_n - k];
    }

    for (k = b_n + 1; k <= n2; k++) {
      costab->data[k] = -costab1q->data[n2 - k];
      sintab->data[k] = costab1q->data[k - b_n];
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
    if (y->size[0] > 1) {
      e = 1.0 / (real_T)y->size[0];
      iDelta2 = y->size[0];
      emxEnsureCapacity((emxArray__common *)y, iDelta2, (int32_T)sizeof(creal_T));
      nd2 = y->size[0] - 1;
      for (iDelta2 = 0; iDelta2 <= nd2; iDelta2++) {
        y->data[iDelta2].re *= e;
        y->data[iDelta2].im *= e;
      }
    }
  }
}

static void eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y)
{
  int32_T n;
  int32_T k;
  int32_T i;
  int32_T j;
  n = x->size[1];
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x->data[i - 1]) {
      k++;
    }
  }

  j = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = k;
  emxEnsureCapacity((emxArray__common *)y, j, (int32_T)sizeof(int32_T));
  j = 0;
  for (i = 1; i <= n; i++) {
    if (x->data[i - 1]) {
      y->data[j] = i;
      j++;
    }
  }
}

void mf_coder(const creal_T rx[16000], real_T tauhat, real_T rolloff, real_T q,
              real_T T, const real_T P[30], const real_T D[3960], real_T RrcNsym,
              emxArray_creal_T *mfout)
{
  real_T Ts;
  int32_T ixstart;
  real_T mtmp;
  int32_T nm1d2;
  boolean_T exitg4;
  real_T ndbl;
  boolean_T exitg3;
  real_T mini;
  boolean_T exitg2;
  boolean_T exitg1;
  real_T maxi;
  real_T y;
  real_T b_y;
  int32_T n;
  real_T apnd;
  real_T cdiff;
  real_T absa;
  real_T absb;
  emxArray_real_T *ns;
  int32_T i7;
  int32_T k;
  emxArray_boolean_T *idx1;
  emxArray_real_T *r10;
  emxArray_real_T *c_y;
  emxArray_real_T *r11;
  emxArray_boolean_T *idx2;
  emxArray_boolean_T *idx3;
  uint32_T varargin_1[2];
  emxArray_real_T *lvec;
  emxArray_int32_T *r12;
  emxArray_int32_T *r13;
  emxArray_real_T *x;
  emxArray_real_T *b_x;
  emxArray_real_T *r14;
  emxArray_real_T *r15;
  emxArray_real_T *z;
  emxArray_real_T *r16;
  static creal_T zvec[16000];
  emxArray_real_T *b_lvec;
  emxArray_creal_T *h;
  emxArray_creal_T *c_x;
  emxArray_creal_T *b_y1;
  emxArray_creal_T *c_y1;
  emxArray_creal_T *d_x;
  emxArray_int32_T *r17;
  emxArray_creal_T *b_h;
  int32_T iv4[2];
  emxArray_int32_T r18;

  /* Run a matched filter on the received samples, returns recieved */
  /* symbols. */
  /*      */
  /*     mfout = mf(rx, tauhat, acqparams) */
  /*  */
  /*  INPUTS: */
  /*    rx              - sampled received signal */
  /*    tauhat        - time offset estimate */
  /*    acqparam        - struct holding aquisition parameters */
  /*        .T              - symbol period */
  /*        .L              - number of symbols (assumed consequtive) */
  /*        .OS             - oversampling factor */
  /*        .pshape         - function representing transmission pulse shape */
  /*  */
  /*  OUTPUTS: */
  /*    mfout        - matched filtered signal at symbol rate */
  /*  (c) 2012 ITR-UniSA */
  /*  Authors: Robby McKilliam */
  /*  Created: 4 April 2012 */
  /* g = acqparam.pshape; */
  /* q = acqparam.OS; */
  /* T = acqparam.T; */
  Ts = mrdivide(T, q);
  ixstart = 1;
  mtmp = P[0];
  if (rtIsNaN(P[0])) {
    nm1d2 = 2;
    exitg4 = 0U;
    while ((exitg4 == 0U) && (nm1d2 < 31)) {
      ixstart = nm1d2;
      if (!rtIsNaN(P[nm1d2 - 1])) {
        mtmp = P[nm1d2 - 1];
        exitg4 = 1U;
      } else {
        nm1d2++;
      }
    }
  }

  if (ixstart < 30) {
    while (ixstart + 1 < 31) {
      if (P[ixstart] < mtmp) {
        mtmp = P[ixstart];
      }

      ixstart++;
    }
  }

  ixstart = 1;
  ndbl = D[0];
  if (rtIsNaN(D[0])) {
    nm1d2 = 2;
    exitg3 = 0U;
    while ((exitg3 == 0U) && (nm1d2 < 3961)) {
      ixstart = nm1d2;
      if (!rtIsNaN(D[nm1d2 - 1])) {
        ndbl = D[nm1d2 - 1];
        exitg3 = 1U;
      } else {
        nm1d2++;
      }
    }
  }

  if (ixstart < 3960) {
    while (ixstart + 1 < 3961) {
      if (D[ixstart] < ndbl) {
        ndbl = D[ixstart];
      }

      ixstart++;
    }
  }

  mini = (mtmp <= ndbl) || rtIsNaN(ndbl) ? mtmp : ndbl;
  ixstart = 1;
  mtmp = P[0];
  if (rtIsNaN(P[0])) {
    nm1d2 = 2;
    exitg2 = 0U;
    while ((exitg2 == 0U) && (nm1d2 < 31)) {
      ixstart = nm1d2;
      if (!rtIsNaN(P[nm1d2 - 1])) {
        mtmp = P[nm1d2 - 1];
        exitg2 = 1U;
      } else {
        nm1d2++;
      }
    }
  }

  if (ixstart < 30) {
    while (ixstart + 1 < 31) {
      if (P[ixstart] > mtmp) {
        mtmp = P[ixstart];
      }

      ixstart++;
    }
  }

  ixstart = 1;
  ndbl = D[0];
  if (rtIsNaN(D[0])) {
    nm1d2 = 2;
    exitg1 = 0U;
    while ((exitg1 == 0U) && (nm1d2 < 3961)) {
      ixstart = nm1d2;
      if (!rtIsNaN(D[nm1d2 - 1])) {
        ndbl = D[nm1d2 - 1];
        exitg1 = 1U;
      } else {
        nm1d2++;
      }
    }
  }

  if (ixstart < 3960) {
    while (ixstart + 1 < 3961) {
      if (D[ixstart] > ndbl) {
        ndbl = D[ixstart];
      }

      ixstart++;
    }
  }

  maxi = (mtmp >= ndbl) || rtIsNaN(ndbl) ? mtmp : ndbl;
  y = q * mini;
  b_y = q * maxi;
  mtmp = q * mini - 16000.0;
  ndbl = q * maxi - 1.0;
  if (rtIsNaN(mtmp) || rtIsNaN(ndbl)) {
    n = 1;
    mtmp = rtNaN;
    apnd = b_y - 1.0;
  } else if (b_y - 1.0 < y - 16000.0) {
    n = 0;
    mtmp = y - 16000.0;
    apnd = b_y - 1.0;
  } else if (rtIsInf(mtmp) || rtIsInf(ndbl)) {
    n = 1;
    mtmp = rtNaN;
    apnd = b_y - 1.0;
  } else {
    mtmp = y - 16000.0;
    ndbl = floor(((b_y - 1.0) - (y - 16000.0)) + 0.5);
    apnd = (y - 16000.0) + ndbl;
    cdiff = apnd - (b_y - 1.0);
    absa = fabs(y - 16000.0);
    absb = fabs(b_y - 1.0);
    if (absa > absb) {
      absb = absa;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
      ndbl++;
      apnd = b_y - 1.0;
    } else if (cdiff > 0.0) {
      apnd = (y - 16000.0) + (ndbl - 1.0);
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int32_T)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&ns, 2);
  i7 = ns->size[0] * ns->size[1];
  ns->size[0] = 1;
  ns->size[1] = n;
  emxEnsureCapacity((emxArray__common *)ns, i7, (int32_T)sizeof(real_T));
  if (n > 0) {
    ns->data[0] = mtmp;
    if (n > 1) {
      ns->data[n - 1] = apnd;
      ixstart = n - 1;
      nm1d2 = ixstart / 2;
      for (k = 1; k <= nm1d2 - 1; k++) {
        ns->data[k] = mtmp + (real_T)k;
        ns->data[(n - k) - 1] = apnd - (real_T)k;
      }

      if (nm1d2 << 1 == ixstart) {
        ns->data[nm1d2] = (mtmp + apnd) / 2.0;
      } else {
        ns->data[nm1d2] = mtmp + (real_T)nm1d2;
        ns->data[nm1d2 + 1] = apnd - (real_T)nm1d2;
      }
    }
  }

  i7 = ns->size[0] * ns->size[1];
  ns->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)ns, i7, (int32_T)sizeof(real_T));
  ixstart = ns->size[0];
  nm1d2 = ns->size[1];
  nm1d2 = ixstart * nm1d2 - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    ns->data[i7] = -ns->data[i7] * Ts - tauhat;
  }

  emxInit_boolean_T(&idx1, 2);
  emxInit_real_T(&r10, 2);

  /* lvec = conj( g(-ns*Ts/p - tauhat) ); */
  /*  root raised cosine pulse */
  /*  */
  /*  INPUT: */
  /*    t       - vector of time instants at which pulse is evaluated */
  /*    rolloff - rolloff factor */
  /*    T       - symbol period */
  /*    OS      - oversampling factor */
  /*  sample period */
  apnd = 1.0 / (2.0 * T);

  /*  cut off frequency */
  cdiff = 1.0 / (T / q);

  /*  sampling frequency */
  /*  implementation based on MATLAB's firrcos */
  /*  use logical indexing: */
  b_abs(ns, r10);
  i7 = idx1->size[0] * idx1->size[1];
  idx1->size[0] = 1;
  idx1->size[1] = r10->size[1];
  emxEnsureCapacity((emxArray__common *)idx1, i7, (int32_T)sizeof(boolean_T));
  nm1d2 = r10->size[0] * r10->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    idx1->data[i7] = (r10->data[i7] < 4.4408920985006262E-16);
  }

  emxInit_real_T(&c_y, 2);

  /*  idx2 = (abs(abs(8*rolloff*fc*t) - 1.0) < sqrt(eps)); */
  y = 8.0 * rolloff * apnd;
  i7 = c_y->size[0] * c_y->size[1];
  c_y->size[0] = 1;
  c_y->size[1] = ns->size[1];
  emxEnsureCapacity((emxArray__common *)c_y, i7, (int32_T)sizeof(real_T));
  nm1d2 = ns->size[0] * ns->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    c_y->data[i7] = y * ns->data[i7];
  }

  emxInit_real_T(&r11, 2);
  b_abs(c_y, r10);
  i7 = r11->size[0] * r11->size[1];
  r11->size[0] = 1;
  r11->size[1] = r10->size[1];
  emxEnsureCapacity((emxArray__common *)r11, i7, (int32_T)sizeof(real_T));
  emxFree_real_T(&c_y);
  nm1d2 = r10->size[0] * r10->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    r11->data[i7] = r10->data[i7] - 1.0;
  }

  emxInit_boolean_T(&idx2, 2);
  b_abs(r11, r10);
  i7 = idx2->size[0] * idx2->size[1];
  idx2->size[0] = 1;
  idx2->size[1] = r10->size[1];
  emxEnsureCapacity((emxArray__common *)idx2, i7, (int32_T)sizeof(boolean_T));
  emxFree_real_T(&r11);
  nm1d2 = r10->size[0] * r10->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    idx2->data[i7] = (r10->data[i7] < 4.4408920985006262E-16);
  }

  emxInit_boolean_T(&idx3, 2);
  i7 = idx3->size[0] * idx3->size[1];
  idx3->size[0] = 1;
  idx3->size[1] = idx1->size[1];
  emxEnsureCapacity((emxArray__common *)idx3, i7, (int32_T)sizeof(boolean_T));
  nm1d2 = idx1->size[0] * idx1->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    idx3->data[i7] = ((!idx1->data[i7]) && (!idx2->data[i7]));
  }

  for (i7 = 0; i7 < 2; i7++) {
    varargin_1[i7] = (uint32_T)ns->size[i7];
  }

  emxInit_real_T(&lvec, 2);
  i7 = lvec->size[0] * lvec->size[1];
  lvec->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)lvec, i7, (int32_T)sizeof(real_T));
  i7 = lvec->size[0] * lvec->size[1];
  lvec->size[1] = (int32_T)varargin_1[1];
  emxEnsureCapacity((emxArray__common *)lvec, i7, (int32_T)sizeof(real_T));
  nm1d2 = (int32_T)varargin_1[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    lvec->data[i7] = 0.0;
  }

  b_emxInit_int32_T(&r12, 2);
  b_emxInit_int32_T(&r13, 2);
  eml_li_find(idx1, r13);
  i7 = r12->size[0] * r12->size[1];
  r12->size[0] = 1;
  r12->size[1] = r13->size[1];
  emxEnsureCapacity((emxArray__common *)r12, i7, (int32_T)sizeof(int32_T));
  nm1d2 = r13->size[0] * r13->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    r12->data[i7] = r13->data[i7];
  }

  mtmp = -sqrt(2.0 * apnd) / (3.1415926535897931 * cdiff);
  nm1d2 = r12->size[0] * r12->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    lvec->data[r12->data[i7] - 1] = mtmp * (3.1415926535897931 * (rolloff - 1.0)
      - 4.0 * rolloff);
  }

  eml_li_find(idx2, r13);
  i7 = r12->size[0] * r12->size[1];
  r12->size[0] = 1;
  r12->size[1] = r13->size[1];
  emxEnsureCapacity((emxArray__common *)r12, i7, (int32_T)sizeof(int32_T));
  emxFree_boolean_T(&idx2);
  nm1d2 = r13->size[0] * r13->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    r12->data[i7] = r13->data[i7];
  }

  y = sqrt(2.0 * apnd) / (6.2831853071795862 * cdiff) * ((3.1415926535897931 *
    (rolloff + 1.0) * sin(3.1415926535897931 * (rolloff + 1.0) / (4.0 * rolloff))
    - 4.0 * rolloff * sin(3.1415926535897931 * (rolloff - 1.0) / (4.0 * rolloff)))
    + 3.1415926535897931 * (rolloff - 1.0) * cos(3.1415926535897931 * (rolloff -
    1.0) / (4.0 * rolloff)));
  nm1d2 = r12->size[0] * r12->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    lvec->data[r12->data[i7] - 1] = y;
  }

  emxInit_real_T(&x, 2);
  eml_li_find(idx3, r13);
  i7 = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = r13->size[1];
  emxEnsureCapacity((emxArray__common *)x, i7, (int32_T)sizeof(real_T));
  ndbl = (1.0 + rolloff) * 2.0 * 3.1415926535897931 * apnd;
  nm1d2 = r13->size[0] * r13->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    x->data[i7] = ndbl * ns->data[r13->data[i7] - 1];
  }

  i7 = r10->size[0] * r10->size[1];
  r10->size[0] = 1;
  r10->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)r10, i7, (int32_T)sizeof(real_T));
  nm1d2 = x->size[0] * x->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    r10->data[i7] = x->data[i7];
  }

  for (k = 0; k <= x->size[1] - 1; k++) {
    r10->data[(int32_T)(1.0 + (real_T)k) - 1] = cos(r10->data[(int32_T)(1.0 +
      (real_T)k) - 1]);
  }

  eml_li_find(idx3, r13);
  i7 = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = r13->size[1];
  emxEnsureCapacity((emxArray__common *)x, i7, (int32_T)sizeof(real_T));
  ndbl = (1.0 - rolloff) * 2.0 * 3.1415926535897931 * apnd;
  nm1d2 = r13->size[0] * r13->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    x->data[i7] = ndbl * ns->data[r13->data[i7] - 1];
  }

  emxInit_real_T(&b_x, 2);
  i7 = b_x->size[0] * b_x->size[1];
  b_x->size[0] = 1;
  b_x->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b_x, i7, (int32_T)sizeof(real_T));
  nm1d2 = x->size[0] * x->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    b_x->data[i7] = x->data[i7];
  }

  for (k = 0; k <= x->size[1] - 1; k++) {
    b_x->data[(int32_T)(1.0 + (real_T)k) - 1] = sin(b_x->data[(int32_T)(1.0 +
      (real_T)k) - 1]);
  }

  emxInit_real_T(&r14, 2);
  eml_li_find(idx3, r13);
  i7 = r14->size[0] * r14->size[1];
  r14->size[0] = 1;
  r14->size[1] = r13->size[1];
  emxEnsureCapacity((emxArray__common *)r14, i7, (int32_T)sizeof(real_T));
  ndbl = 8.0 * rolloff * apnd;
  nm1d2 = r13->size[0] * r13->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    r14->data[i7] = ndbl * ns->data[r13->data[i7] - 1];
  }

  rdivide(b_x, r14, x);
  eml_li_find(idx3, r13);
  i7 = r12->size[0] * r12->size[1];
  r12->size[0] = 1;
  r12->size[1] = r13->size[1];
  emxEnsureCapacity((emxArray__common *)r12, i7, (int32_T)sizeof(int32_T));
  emxFree_real_T(&r14);
  nm1d2 = r13->size[0] * r13->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    r12->data[i7] = r13->data[i7];
  }

  emxInit_real_T(&r15, 2);
  mtmp = -4.0 * rolloff / cdiff;
  eml_li_find(idx3, r13);
  i7 = r15->size[0] * r15->size[1];
  r15->size[0] = 1;
  r15->size[1] = r13->size[1];
  emxEnsureCapacity((emxArray__common *)r15, i7, (int32_T)sizeof(real_T));
  ndbl = 8.0 * rolloff * apnd;
  emxFree_boolean_T(&idx3);
  nm1d2 = r13->size[0] * r13->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    r15->data[i7] = ndbl * ns->data[r13->data[i7] - 1];
  }

  emxFree_int32_T(&r13);
  emxInit_real_T(&z, 2);
  power(r15, b_x);
  i7 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = r10->size[1];
  emxEnsureCapacity((emxArray__common *)z, i7, (int32_T)sizeof(real_T));
  emxFree_real_T(&r15);
  nm1d2 = r10->size[0] * r10->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    z->data[i7] = mtmp * (r10->data[i7] + x->data[i7]);
  }

  emxInit_real_T(&r16, 2);
  i7 = r16->size[0] * r16->size[1];
  r16->size[0] = 1;
  r16->size[1] = b_x->size[1];
  emxEnsureCapacity((emxArray__common *)r16, i7, (int32_T)sizeof(real_T));
  ndbl = 3.1415926535897931 * sqrt(1.0 / (2.0 * apnd));
  nm1d2 = b_x->size[0] * b_x->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    r16->data[i7] = ndbl * (b_x->data[i7] - 1.0);
  }

  emxFree_real_T(&b_x);
  rdivide(z, r16, r10);
  emxFree_real_T(&r16);
  emxFree_real_T(&z);
  nm1d2 = r10->size[0] * r10->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    lvec->data[r12->data[i7] - 1] = r10->data[i7];
  }

  emxFree_int32_T(&r12);

  /*  scale pulse */
  mtmp = sqrt(2.0 * apnd);
  i7 = lvec->size[0] * lvec->size[1];
  lvec->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)lvec, i7, (int32_T)sizeof(real_T));
  ixstart = lvec->size[0];
  nm1d2 = lvec->size[1];
  nm1d2 = ixstart * nm1d2 - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    lvec->data[i7] *= mtmp;
  }

  y = RrcNsym * q * Ts;
  b_abs(ns, r10);
  i7 = idx1->size[0] * idx1->size[1];
  idx1->size[0] = 1;
  idx1->size[1] = r10->size[1];
  emxEnsureCapacity((emxArray__common *)idx1, i7, (int32_T)sizeof(boolean_T));
  emxFree_real_T(&ns);
  nm1d2 = r10->size[0] * r10->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    idx1->data[i7] = (r10->data[i7] <= y);
  }

  emxFree_real_T(&r10);
  i7 = lvec->size[0] * lvec->size[1];
  lvec->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)lvec, i7, (int32_T)sizeof(real_T));
  ixstart = lvec->size[0];
  nm1d2 = lvec->size[1];
  nm1d2 = ixstart * nm1d2 - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    lvec->data[i7] *= (real_T)idx1->data[i7];
  }

  emxFree_boolean_T(&idx1);
  for (i7 = 0; i7 < 16000; i7++) {
    zvec[i7].re = 0.0;
    zvec[i7].im = 0.0;
    zvec[i7] = rx[i7];
  }

  /* zerofilled r */
  /* compute a convolution using the fft.  The output is identical to: */
  /* conv(a,b,`valid')  when the length(a) >= length(b),  and */
  /* conv(b,a,`valid') when length(a) =< length(b). */
  /*  */
  /*  (c) 2012 ITR-UniSA */
  /*  Authors: Robby McKilliam */
  /*  Created: 3 April 2012 */
  varargin_1[0] = (uint32_T)lvec->size[1];
  varargin_1[1] = 16000U;
  n = (int32_T)varargin_1[0];
  if (16000 < (int32_T)varargin_1[0]) {
    n = 16000;
  }

  b_emxInit_real_T(&b_lvec, 1);
  i7 = b_lvec->size[0];
  b_lvec->size[0] = lvec->size[1];
  emxEnsureCapacity((emxArray__common *)b_lvec, i7, (int32_T)sizeof(real_T));
  nm1d2 = lvec->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    b_lvec->data[i7] = lvec->data[i7];
  }

  emxInit_creal_T(&h, 2);
  emxInit_creal_T(&c_x, 2);
  b_emxInit_creal_T(&b_y1, 1);
  b_emxInit_creal_T(&c_y1, 1);
  b_eml_fft(b_lvec, (uint32_T)lvec->size[1] + 15999U, b_y1);
  c_eml_fft(zvec, (uint32_T)lvec->size[1] + 15999U, c_y1);
  ipermute(b_y1, c_x);
  ipermute(c_y1, h);
  i7 = c_x->size[0] * c_x->size[1];
  c_x->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)c_x, i7, (int32_T)sizeof(creal_T));
  ixstart = c_x->size[0];
  nm1d2 = c_x->size[1];
  emxFree_real_T(&b_lvec);
  emxFree_creal_T(&c_y1);
  emxFree_real_T(&lvec);
  nm1d2 = ixstart * nm1d2 - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    mtmp = c_x->data[i7].re;
    ndbl = c_x->data[i7].im;
    apnd = h->data[i7].re;
    cdiff = h->data[i7].im;
    c_x->data[i7].re = mtmp * apnd - ndbl * cdiff;
    c_x->data[i7].im = mtmp * cdiff + ndbl * apnd;
  }

  b_emxInit_creal_T(&d_x, 1);
  i7 = d_x->size[0];
  d_x->size[0] = c_x->size[1];
  emxEnsureCapacity((emxArray__common *)d_x, i7, (int32_T)sizeof(creal_T));
  nm1d2 = c_x->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    d_x->data[i7] = c_x->data[i7];
  }

  i7 = c_x->size[1];
  if (i7 < 0) {
    i7 = 0;
  }

  d_eml_fft(d_x, (uint32_T)i7, b_y1);
  ipermute(b_y1, h);
  ndbl = (real_T)(h->size[1] - n) + 1.0;
  emxFree_creal_T(&d_x);
  emxFree_creal_T(&b_y1);
  emxFree_creal_T(&c_x);
  if ((real_T)n > ndbl) {
    n = 1;
    i7 = 0;
  } else {
    i7 = (int32_T)ndbl;
  }

  emxInit_int32_T(&r17, 1);
  k = r17->size[0];
  r17->size[0] = (i7 - n) + 1;
  emxEnsureCapacity((emxArray__common *)r17, k, (int32_T)sizeof(int32_T));
  nm1d2 = i7 - n;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    r17->data[i7] = n + i7;
  }

  emxInit_creal_T(&b_h, 2);
  iv4[0] = 1;
  iv4[1] = r17->size[0];
  i7 = b_h->size[0] * b_h->size[1];
  b_h->size[0] = iv4[0];
  b_h->size[1] = iv4[1];
  emxEnsureCapacity((emxArray__common *)b_h, i7, (int32_T)sizeof(creal_T));
  nm1d2 = iv4[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    ixstart = iv4[0] - 1;
    for (k = 0; k <= ixstart; k++) {
      r18 = *r17;
      r18.size = (int32_T *)&iv4;
      r18.numDimensions = 1;
      b_h->data[k + b_h->size[0] * i7] = h->data[r18.data[k + r18.size[0] * i7]
        - 1];
    }
  }

  emxFree_int32_T(&r17);
  i7 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = b_h->size[1];
  emxEnsureCapacity((emxArray__common *)h, i7, (int32_T)sizeof(creal_T));
  nm1d2 = b_h->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    h->data[h->size[0] * i7] = b_h->data[b_h->size[0] * i7];
  }

  emxFree_creal_T(&b_h);

  /* plot(real(h)); */
  /* length(h) */
  /* length(q*(mini:maxi)) */
  if (rtIsNaN(mini) || rtIsNaN(maxi)) {
    n = 1;
    mtmp = rtNaN;
    apnd = maxi;
  } else if (maxi < mini) {
    n = 0;
    mtmp = mini;
    apnd = maxi;
  } else if (rtIsInf(mini) || rtIsInf(maxi)) {
    n = 1;
    mtmp = rtNaN;
    apnd = maxi;
  } else {
    mtmp = mini;
    ndbl = floor((maxi - mini) + 0.5);
    apnd = mini + ndbl;
    cdiff = apnd - maxi;
    absa = fabs(mini);
    absb = fabs(maxi);
    if (absa > absb) {
      absb = absa;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
      ndbl++;
      apnd = maxi;
    } else if (cdiff > 0.0) {
      apnd = mini + (ndbl - 1.0);
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int32_T)ndbl;
    } else {
      n = 0;
    }
  }

  i7 = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = n;
  emxEnsureCapacity((emxArray__common *)x, i7, (int32_T)sizeof(real_T));
  if (n > 0) {
    x->data[0] = mtmp;
    if (n > 1) {
      x->data[n - 1] = apnd;
      ixstart = n - 1;
      nm1d2 = ixstart / 2;
      for (k = 1; k <= nm1d2 - 1; k++) {
        x->data[k] = mtmp + (real_T)k;
        x->data[(n - k) - 1] = apnd - (real_T)k;
      }

      if (nm1d2 << 1 == ixstart) {
        x->data[nm1d2] = (mtmp + apnd) / 2.0;
      } else {
        x->data[nm1d2] = mtmp + (real_T)nm1d2;
        x->data[nm1d2 + 1] = apnd - (real_T)nm1d2;
      }
    }
  }

  y = q * mini;
  i7 = mfout->size[0] * mfout->size[1];
  mfout->size[0] = 1;
  mfout->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)mfout, i7, (int32_T)sizeof(creal_T));
  nm1d2 = x->size[0] * x->size[1] - 1;
  for (i7 = 0; i7 <= nm1d2; i7++) {
    mfout->data[i7].re = Ts * h->data[(int32_T)((q * x->data[i7] - y) + 1.0) - 1]
      .re;
    mfout->data[i7].im = Ts * h->data[(int32_T)((q * x->data[i7] - y) + 1.0) - 1]
      .im;
  }

  emxFree_real_T(&x);
  emxFree_creal_T(&h);
}

/* End of code generation (mf_coder.cpp) */
