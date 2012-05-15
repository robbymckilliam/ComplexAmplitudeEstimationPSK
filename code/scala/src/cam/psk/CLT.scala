/*
 * Classes for computing the central limit theorem of the least
 * squares phase and amplitude estimator
 * 
 * @author Robby McKilliam
 */

package cam.psk

import pubsim.distributions.complex.ComplexRandomVariable
import numbers.finite.integration.RealIntegral

trait CLT {
  /** The variance of the phase and amplitude as a tuple (phase, amp) */
  def variance : (Double, Double)
  /** The covariance of the phase and amplitude */
  def covariance  : Double
}

/** Some convenience functions are defined here */
abstract class AbstractCLT extends CLT {
  def sin(x : Double) = scala.math.sin(x)
  def cos(x : Double) = scala.math.cos(x)
  def sqr(x : Double) = x*x
  def sqrt(x : Double) = scala.math.sqrt(x)
  val pi = scala.math.Pi
}

/** 
 * Computes the central limit theorem for the coherent version of the least
 * square estimator.  The proportion of pilots symbols is p which must be in
 * the interval [0,1].  The proportion of data symbols is correspondinly 
 * d = 1 - p.
 * 
 * X is a complex random variable that is assumed to be circularly symmetric
 * and satisfy all of the properties required by Theorem 1.
 * 
 */
class CoherentCLT(M : Int, X : ComplexRandomVariable, p : Double) extends AbstractCLT {
 
  val Z = X.magnitudeMarginal
  
  protected def fZ(z : Double) = Z.pdf(z)
  
  protected def f(r : Double, phi : Double) : Double = {
    val z = sqrt(r*r - 2*r*cos(phi) + 1)
    return r*fZ(z)/z
  }
  
  protected def g(phi : Double) = RealIntegral.trapezoidal( r => f(r,phi), 0, 20*Z.getVariance, 1000) 
  
  protected def fracpart(x : Double) = x - 2*pi/M*scala.math.round(M*x/2/pi)
  
  val A1 = RealIntegral.trapezoidal( x => sqr(sin(x))*g(x), 0, 2*pi, 1000)
  val A2 = RealIntegral.trapezoidal( x => sqr(sin(fracpart(x)))*g(x), 0, 2*pi, 1000)
  val H =  2*sin(pi/M) * (0 to M-1).map(k => g((2*k+1)*pi/M)).foldLeft(0.0){ (s : Double ,v : Double) => s + v }
  val d = 1 - p
  
  override def variance : (Double, Double) = ( (A1 + A2)/sqr(p + d*H), 0.0 )
  
  override def covariance  : Double = 0.0
  
}


