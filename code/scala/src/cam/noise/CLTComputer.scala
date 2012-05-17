/*
 * Classes that wraps noise distributions with functions relavant to computing the asymptotic
 * variance the least squares estimator.
 * @author Robby McKilliam
 */

package cam.noise

import numbers.finite.integration.RealIntegral
import pubsim.distributions.complex.ComplexRandomVariable
import pubsim.distributions.Chi

/** Function names come directly from the paper */
trait CLTComputer {
  
  /** The marginal pdf of the magnitude of the complex random variable */
  def fZ(z : Double) : Double
  /** The joint pdf of magnitude and phase of the the translated complex random variable*/
  def f(r : Double, phi : Double) : Double = {
    val z = sqrt(r*r - 2*r*cos(phi) + 1)
    return r*fZ(z)/z/2/pi
  }
  /** The marginal pdf of the phase of the translated complex random variable */
  def g(phi : Double) : Double
  def g2(phi : Double) : Double
  def h1(x : Double) : Double
  def h2(x : Double) : Double
  def G(x : Double) : Double
  
  val A1 : Double
  val A2 : Double
  val H : Double
  
  /** Triple returns phase variance, amplitude variance and phase amplitude covariance */
  def variance : (Double, Double, Double)
  
  // Some convenience functions are defined here
  def sin(x : Double) = scala.math.sin(x)
  def cos(x : Double) = scala.math.cos(x)
  def sqr(x : Double) = x*x
  def sqrt(x : Double) = scala.math.sqrt(x)
  val pi = scala.math.Pi
  
}

abstract class AbstractCLTComputer(M : Int, p : Double) extends CLTComputer {
  protected val d = 1 - p
  def fracpart(x : Double) = x - 2*pi/M*scala.math.round(M*x/2/pi)
  
  override lazy val A1 = RealIntegral.trapezoidal( x => sqr(sin(x))*g2(x), -pi, pi, 1000)
  override lazy val A2 = RealIntegral.trapezoidal( x => sqr(sin(fracpart(x)))*g2(x), -pi, pi, 1000)
  override def h2(x : Double) = RealIntegral.trapezoidal( phi => cos(fracpart(x + phi))*g(phi), -pi, pi, 1000)
  override def h1(x : Double) = RealIntegral.trapezoidal( phi => cos(x + phi)*g(phi), -pi, pi, 1000)
  override lazy val H =  h2(0) - 2*sin(pi/M) * (0 to M-1).map(k => g((2*k+1)*pi/M)).foldLeft(0.0){ (s : Double ,v : Double) => s + v }
  override def G(x : Double) = p*h1(x) + d*h2(x)
  
  override def variance : (Double, Double, Double) = {
    ((p*A1+d*A2)/sqr(p + H*d), 0.0, 0.0)
  }
  
}

/** 
 * Takes a ComplexRandomVariable and computes all the necessary 
 * values by numerical integration.  Might not be numerically accurate and is fairly slow.
 * THIS CURRENTLY ASSUMES THAT the true amplitude rho0 = 1
 */
class GeneralCLTComputer(M : Int, p : Double, X : ComplexRandomVariable) extends AbstractCLTComputer(M,p) {
  
  //get the class representing the marginal magnitude random variable of X
  protected val Z = X.magnitudeMarginal
  
  override def fZ(z : Double) = Z.pdf(z)
  
  //Compute g by numerical integration
  override def g(phi : Double) = RealIntegral.trapezoidal( r => r*f(r,phi), 0, 30*sqrt(Z.getVariance), 1000) 
  
  //Compute g by numerical integration
  override def g2(phi : Double) = RealIntegral.trapezoidal( r => r*r*f(r,phi), 0, 30*sqrt(Z.getVariance), 1000) 
  
  /** The marginal pdf of the phase */
  def f(phi : Double) : Double = RealIntegral.trapezoidal( r => f(r,phi), 0, 30*sqrt(Z.getVariance), 1000) 
  
}

/** 
 *CLT for circularly symmetric complex Gaussian noise.  Makes use of specific formula for
 *this case.  This is fast and should be numerically accurate.
*/
class GaussianCLT(M : Int, p : Double, sigma : Double, rho0 : Double) extends AbstractCLTComputer(M,p) {
  
  /** Signal amplitude defaults to 1 */
  def this(M : Int, p: Double, sigma : Double) = this(M,p,sigma,1.0)
  
  protected val k = rho0/sigma
  protected val Z = new Chi.Chi2(1/k/k)   //marginal distribution of the magnitude
  
  override def fZ(z : Double) = Z.pdf(z)
  
  //Compute g by numerical integration
  override def g(phi : Double) = RealIntegral.trapezoidal( r => r*f(r,phi), 0, 30*sqrt(Z.getVariance), 1000) 
  
  //Compute g by numerical integration
  override def g2(phi : Double) = RealIntegral.trapezoidal( r => r*r*f(r,phi), 0, 30*sqrt(Z.getVariance), 1000) 
  
  /** The marginal pdf of the phase */
  def f(phi : Double) : Double = RealIntegral.trapezoidal( r => f(r,phi), 0, 30*sqrt(Z.getVariance), 1000) 
  
}