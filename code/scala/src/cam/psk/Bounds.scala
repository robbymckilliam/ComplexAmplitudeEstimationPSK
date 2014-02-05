/*
 * Classes for computing estimator bounds for noncoherent PSK carrier phase
 * estimation
 * 
 * @author Robby McKilliam
 */

package cam.psk

import numbers.finite.integration.RealIntegral
import scala.math.exp
import scala.math.Pi
import scala.math.sqrt
import scala.math.sinh
import scala.math.cosh
import scala.math.tanh

trait CramerRaoBound {
  /** Returns a bound on the estimator variance given noise variance v and number of symbols L*/
  def variance(v : Double, L : Int) : Double
}

/** 
 * This is the bound derived by Bill Cowley for BPSK
 */
class BPSKCRB extends CramerRaoBound {
  
  final def variance(v : Double, L : Int) : Double = {  
    val Fb = RealIntegral.trapezoidal( n => f(v,n), -30.0*v, 30.0*v, 100000)
    return v/L/Fb
  }
  
  private def f(v : Double, n : Double) : Double = {
    val tan1 = tanh((1+n)/v)
    val exp1 = exp(-n*n/2/v)/sqrt(2*Pi*v)
    return tan1*tan1*exp1
  }
}

/** 
 * This is the bound derived by Bill Cowley for QPSK
 */
class QPSKCRB extends CramerRaoBound {
  
  final def variance(v : Double, L : Int) : Double = {  
    def Fqi(nq : Double) = RealIntegral.trapezoidal( ni => f(v,ni, nq), -20.0*sqrt(v), 20.0*sqrt(v), 1000)
    val Fq = RealIntegral.trapezoidal( nq => Fqi(nq), -20.0*sqrt(v), 20.0*sqrt(v), 1000)
    return v/L/Fq
  }
  
  private def f(v : Double, ni : Double, nq : Double) : Double = {
    val sinhi = sinh((1+ni)/v)*nq
    val sinhq = sinh(nq/v)*(1+ni)
    val cosh1 = cosh((1+ni)/v) + cosh(nq/v)
    val exp1 = scala.math.exp(-(ni*ni+nq*nq)/2/v)/(2*scala.math.Pi*v)
    val a = (sinhi - sinhq)/cosh1
    return a * a * exp1 / v
  }
}
