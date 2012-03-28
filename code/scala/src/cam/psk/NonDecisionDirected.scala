/*
 * NonDecisionDirected.scala
 * Contains various so called NonDecisionDirected algorithms.  To my knowledge, all of these
 * are based on taking the phase of the recieved symbols to the Mth power, and then doing
 * `something' with the magnitude.  The original idea appears to come from the paper
 * 
 * A. J. Viterbi and A. M. Viterbi "Nonlinear estimation of PSK-modulated carrier phase 
 *  with application to burst digital transmission", IEEE Trans. Inf. Theory
 * 
 * @author Robby McKilliam
 */

package cam.psk

import numbers.finite.Complex
import numbers.finite.PolarComplex

/** 
 * Implement the Viterbi and Viterbi algorithm.
 * Contructor takes:
 * M, as in M-PSK,
 * D, a sequence of indices for the recieved symbols
 * F, a function for modulating the magnitude
 */
class ViterbiViterbi(M : Int, D : Seq[Int], F : Double => Double) 
extends NonCoherentComplexAmplitudeEstimator(M) {
  
  def estimate(y : Seq[Complex]) : Complex = {
    return y.foldLeft(Complex.zero)( (s,v) => s + new PolarComplex(F(v.magnitude),M*v.angle) )
  }
  
}
