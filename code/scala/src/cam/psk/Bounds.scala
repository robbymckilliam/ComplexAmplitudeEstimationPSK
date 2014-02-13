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
import numbers.finite.RectComplex
import scala.math.exp
import scala.math.Pi
import scala.math.sqrt
import scala.math.sinh
import scala.math.cosh
import scala.math.tanh
import scala.math.tanh
import scala.math.Pi

trait CramerRaoBound {
  /** Returns a bound on the variance of the estimator given number of symbols L*/
  def variance(L : Int) : Double
}

/** 
 * This is the bound derived by Bill Cowley for BPSK
 * Constructor takes noise variance
 */
class BPSKCRB(val v : Double)  extends CramerRaoBound {
  
  val Fb = RealIntegral.trapezoidal( n => f(n), -30.0*v, 30.0*v, 100000)
  
  final def variance(L : Int) : Double = {  
    return v/L/Fb
  }
  
  def f(n : Double) : Double = {
    val tan1 = tanh((1+n)/v)
    val exp1 = exp(-n*n/2/v)/sqrt(2*Pi*v)
    return tan1*tan1*exp1
  }
  
}

/** 
 * This is the bound derived by Bill Cowley for QPSK
 * Constructor takes noise variance
 */
class QPSKCRB(val v : Double) extends CramerRaoBound {
  
  def Fqi(nq : Double) = RealIntegral.trapezoidal( ni => f(ni, nq), -20.0*sqrt(v), 20.0*sqrt(v), 1000)
  val Fq = RealIntegral.trapezoidal( nq => Fqi(nq), -20.0*sqrt(v), 20.0*sqrt(v), 1000)
  
  final def variance(L : Int) : Double = {  
    return v/L/Fq
  }
  
  def f(ni : Double, nq : Double) : Double = {
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
 * 
 * Constructor also takes noise variance
 * 
 * The notation here differs from that in the paper above.  There's a description of this in the
 * comments of the tex file.
 * 
 */
abstract class QAMCRB(val q1 : Seq[Complex], val v : Double) extends CramerRaoBound {
  
  //The function Lambda from the paper.  Although it's not made clear in the paper, this is
  //actually a function of phi (the `true' phase). It likely the CRB is invariant to phi, so this 
  //ultimately cause no problem and I've set phi to zero (which is what ultimately is done in
  //the paper too).  I'm not a huge fan of the way it is written in the paper,
  def Lambda(x : Complex) : Double = {    
    val Ax = A(x); val Dx = D(x);
    return -Ax*Ax/v/Dx/Dx + B(x)/Dx/v - C(x)/Dx
  }
  
  def A(x : Complex) = q1.foldLeft(0.0) { (sum, a) => 
    val ax = a*x; val acx = a.conjugate * x;
    sum + s(ax)*ax.imag + s(acx)*acx.imag
  }
  
  def B(x : Complex) = q1.foldLeft(0.0) { (sum, a) => 
    val ax = a*x; val acx = a.conjugate * x;
    sum + c(ax)*ax.imag*ax.imag + c(acx)*acx.imag*acx.imag
  }
  
  def C(x : Complex) = q1.foldLeft(0.0) { (sum, a) => 
    val ax = a*x; val acx = a.conjugate * x;
    sum + s(ax)*ax.real + s(acx)*acx.real
  }
  
  def D(x : Complex) = q1.foldLeft(0.0) { (sum, a) => 
    val ax = a*x; val acx = a.conjugate * x;
    sum + c(ax) + c(acx)
  }
  
  def s(x : Complex) = sinh(x.real/v)
  
  def c(x : Complex) = cosh(x.real/v)

} 

object MPSKCRB {
  def constellation(M: Int) = (0 until M).map( m => PolarComplex(1, m*2*Pi/M))
  def firstquadrant(M: Int) = (0 until M/4).map( m => PolarComplex(1, m*2*Pi/M + Pi/M))
}

class MPSKCRB(val M: Int, override val v: Double) extends QAMCRB(MPSKCRB.firstquadrant(M), v)   {
  
  def Fqi(nq : Double) = RealIntegral.trapezoidal( ni => f(ni, nq), -20.0*sqrt(v), 20.0*sqrt(v), 1000)
  val Fq = RealIntegral.trapezoidal( nq => Fqi(nq), -20.0*sqrt(v), 20.0*sqrt(v), 1000)
  
  final def variance(L : Int) : Double = {  
    return v/L/Fq
  }
  
  def f(ni : Double, nq : Double) : Double = {
    val exp1 = exp(-(ni*ni+nq*nq)/2/v)/(2*scala.math.Pi*v) //Bivariate Guassian pdf
    return Lambda(q1(0)+RectComplex(ni,nq)) * exp1
  }
  
}
