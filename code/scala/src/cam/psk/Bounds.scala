/*
 * Classes for computing estimator bounds for noncoherent PSK carrier phase
 * estimation
 * 
 * @author Robby McKilliam
 */

package cam.psk

import numbers.finite.integration.RealIntegral

trait CramerRaoBound {
  /** Returns a bound on the estimator variance given noise variance v and number of symbols L*/
  def variance(v : Double, L : Int) : Double
}

/** 
 * This is the bound derived by Bill Cowley.
 */
class BPSKCRB extends CramerRaoBound {
  
  final def variance(v : Double, L : Int) : Double = {  
    return RealIntegral.trapezoidal( n => f(v,n), -100.0, 100.0, 1000000) 
  }
  
  private def f(v : Double, n : Double) : Double = {
    val tan = scala.math.tanh((1+n)/v)
    val exp = scala.math.exp(-n*n/2/v)
    return tan*tan*exp/scala.math.sqrt(2*scala.math.Pi*v)
  }
}
