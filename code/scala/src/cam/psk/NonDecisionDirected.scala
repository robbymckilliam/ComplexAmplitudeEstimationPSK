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
class ViterbiViterbi(M : Int, val D : Seq[Int], val F : Double => Double) 
extends NonCoherentComplexAmplitudeEstimator(M) {
  
  val L = D.length
  
  /** The magnitude estimator returned is not accurate! */ 
  def estimate(y : Seq[Complex]) : Complex = {
    val Y = y.foldLeft(Complex.zero)( (s,v) => s + new PolarComplex(F(v.magnitude),M*v.angle) )/L
    return new PolarComplex(Y.magnitude, Y.angle/M)
  }
  
}
