/*
 * dcblockST2DL.cpp
 *
 * Code generation for function 'dcblockST2DL'
 *
 * C source code generated on: Thu Nov  8 13:52:45 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dcblockST2DL.h"
#include "freqest_periodogram_coder.h"
#include "mf_coder.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void dcblockST2DL(const creal_T RxSamp[16000], creal_T RxSampFilt[16000])
{
  creal_T dbuffer[2];
  int32_T j;
  int32_T k;

  /*  DC blocking filter */
  /* ------------------------------------------------------------- */
  /*  dcblock.m  */
  /*  */
  /*  [RxSampFilt] = dcblockST2DL(RxSamp) */
  /*  */
  /*  This digital filter blocks the DC component in the received signal. It is */
  /*  primarily intended to remove the complex exponential at DC that is */
  /*  transmitted in the ST2 downlink waveform for frequency estimation */
  /*  purposes. Once the frequency has been estimated, the discrete tone should */
  /*  be removed in order to avoid bit error rate degradations. The */
  /*  implementation uses an order-one IIR filter. */
  /*   */
  /*  Inputs: */
  /*    RxSamp      - Vector containing samples of received signal */
  /*  */
  /*  Outputs: */
  /*    RxSampFilt  - Vector containing filtered received samples */
  /*  */
  /*  Comments: For details see [J. de Freitas, "The DC Blocking Filter", Jan. */
  /*            2007, available at */
  /*            http://www.mathworks.com.au/matlabcentral/fileexchange/13792-the-dc-blocking-filter]  */
  /*   */
  /* ------------------------------------------------------------- */
  /*  Author: A. Pollok, Institute for Telecommunications Research */
  /*  Project: ASRP */
  /* ------------------------------------------------------------- */
  /*  Copyright 2012  */
  /*  Institute for Telecommunications Research */
  /*  University of South Australia */
  /* ------------------------------------------------------------- */
  /*  ST1 UL: coefficents of order-one IIR filter */
  /*  coeff for cut-on frequency at 12Hz */
  /*  fc = 12 / (OS/(2*Ts)); % normalised cut-on frequency */
  /*  p = dcblock_MatlabCentral(fc); */
  /*  filter frequency corrected samples */
  dbuffer[1].re = 0.0;
  dbuffer[1].im = 0.0;
  for (j = 0; j < 16000; j++) {
    dbuffer[0] = dbuffer[1];
    dbuffer[1].re = 0.0;
    dbuffer[1].im = 0.0;
    for (k = 0; k < 2; k++) {
      dbuffer[k].re += RxSamp[j].re * (1.0 + -2.0 * (real_T)k);
      dbuffer[k].im += RxSamp[j].im * (1.0 + -2.0 * (real_T)k);
    }

    dbuffer[1].re -= dbuffer[0].re * -0.998009283733233;
    dbuffer[1].im -= dbuffer[0].im * -0.998009283733233;
    RxSampFilt[j] = dbuffer[0];
  }
}

/* End of code generation (dcblockST2DL.cpp) */
