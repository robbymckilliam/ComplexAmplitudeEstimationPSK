/*
 * Classes for computing estimator bounds for noncoherent PSK carrier phase
 * estimation
 * 
 * @author Robby McKilliam
 */

package cam.psk

import numbers.finite.integration.RealIntegral
import numbers.finite.Complex
import numbers.finite.PolarComplex
import scala.math.exp
import scala.math.Pi
import scala.math.sqrt
import scala.math.sinh
import scala.math.cosh
import scala.math.tanh
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
    val exp1 = exp(-(ni*ni+nq*nq)/2/v)/(2*scala.math.Pi*v)
    val a = (sinhi - sinhq)/cosh1
    return a * a * exp1 / v
  }
}

/**
 * This is the CRB for QAM, 8PSK etc developed by
 * 
 * Rice, Cowley, Moran, Rice, "Cramér–Rao Lower Bounds for QAM Phase and
 * Frequency Estimation", IEEE TRANSACTIONS ON COMMUNICATIONS, VOL. 49, NO. 9, SEPTEMBER 2001
 *   
 * At construction this requires a set of constellation points from the first quadrant.  Rice et. al. are not
 * specific about what happens to points on the boundary, i.e. on the positive real or imaginary axis.
 * This isn't a problem for appropriately arranged constellations such as 16QAM, 64QAM and 8PSK
 */
class QAMCRB(val q1 : Seq[Complex]) extends CramerRaoBound {
  
   final def variance(v : Double, L : Int) : Double = {  
    return 0.0
  }
  
  //The function Lambda from the paper.  Although it's not made clear in the paper, this is
  //actually a function of phi (the `true' phase). It likely the CRB is invariant to phi, so this 
  //ultimately cause no problem. I'm not a huge fan of the way it is written in the paper though.
  def Lambda(x : Complex, v : Double, phi : Double = 0.0) : Double = {    
    
    //function for computing sums inside all of the terms
    def sumr(f : (Complex, Int) => Double, a : Complex) : Double = {
      val sumr1 = List(-1,1).foldLeft(0.0) { (sum, r) => sum + f(a,r) }
      return exp(-a.mag2/2/v) * sumr1
    }
    
    //denomimator for all terms
    def coshr(a : Complex, r : Int) : Double = {
      val eps = epsilon(phi,x,a,r)
      return cosh(eps.real/v)
    }
    val denom = q1.foldLeft(0.0) { (sum, a) => sum + sumr(coshr,a) }    
    
    //numerator of the first term in the definition of Lambda
    def sinhri(a : Complex, r : Int) : Double = {
      val eps = epsilon(phi,x,a,r)
      return sinh(eps.real/v) * eps.imag
    }
    val num1 = q1.foldLeft(0.0) { (sum, a) => sum + sumr(sinhri,a) }

     //numerator of the first term in the definition of Lambda
    def coshri(a : Complex, r : Int) : Double = {
      val eps = epsilon(phi,x,a,r)
      return cosh(eps.real/v) * eps.imag * eps.imag
    }
    val num2 = q1.foldLeft(0.0) { (sum, a) => sum + sumr(coshri,a) }
    
    return 
    
  }
  
  def epsilon(phi: Double, x : Complex, a : Complex, r : Int) : Complex = {
    val ar = if(r == 1) a.conjugate else if(r == -1) a else throw new RuntimeException("r must be either -1 or 1")
    return ar*x*PolarComplex(1,-phi)
  }

} 
