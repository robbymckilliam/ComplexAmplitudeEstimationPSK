/*
 * Classes that wrap noise distributions with functions relavant to computing the asymptotic
 * variance the least squares estimator.
 * @author Robby McKilliam
 */

package cam.noise

import numbers.finite.integration.RealIntegral.trapezoidal
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
  def m1(x : Double) : Double
  def m2(x : Double) : Double
  def G(x : Double) : Double
  
  val A1 : Double
  val A2 : Double
  val B1 : Double
  val B2 : Double
  val H : Double
  
  /** Triple returns phase variance and ampltiude variance */
  def variance : (Double, Double)
  
  // Some convenience functions are defined here
  final def sin(x : Double) = scala.math.sin(x)
  final def cos(x : Double) = scala.math.cos(x)
  final def sqr(x : Double) = x*x
  final def cub(x : Double) = x*x*x
  final def sqrt(x : Double) = scala.math.sqrt(x)
  final def exp(x : Double) = scala.math.exp(x)
  final def erf(x : Double) = pubsim.Util.erf(x)
  val pi = scala.math.Pi
  
}

abstract class AbstractCLTComputer(val M : Int, val p : Double) extends CLTComputer {
  protected val d = 1 - p
  def fracpart(x : Double) = x - 2*pi/M*scala.math.round(M*x/2/pi)
  
  override lazy val A1 = trapezoidal( x => sqr(sin(x))*g2(x), -pi, pi, 2000) - sqr(m1(0))
  override lazy val A2 = trapezoidal( x => sqr(sin(fracpart(x)))*g2(x), -pi, pi, 2000) - sqr(m2(0))
  override lazy val B1 = trapezoidal( x => sqr(cos(x))*g2(x), -pi, pi, 2000) - 1.0
  override lazy val B2 = trapezoidal( x => sqr(cos(fracpart(x)))*g2(x), -pi, pi, 2000) - sqr(h2(0))
  override def h2(x : Double) = trapezoidal( phi => cos(fracpart(x + phi))*g(phi), -pi, pi, 2000)
  override def h1(x : Double) = trapezoidal( phi => cos(x + phi)*g(phi), -pi, pi, 2000)
  override def m2(x : Double) = trapezoidal( phi => sin(fracpart(x + phi))*g(phi), -pi, pi, 2000)
  override def m1(x : Double) = trapezoidal( phi => sin(x + phi)*g(phi), -pi, pi, 2000)
  override lazy val H =  h2(0) - 2*sin(pi/M) * (0 to M-1).map(k => g((2*k+1)*pi/M)).foldLeft(0.0){ (s : Double ,v : Double) => s + v }
  override def G(x : Double) = p*h1(x) + d*h2(x)
  
  override def variance : (Double, Double) = {
    ((p*A1+d*A2)/sqr(p + H*d), p*B1 + d*B2)
  }
  
}

/** 
 * Takes a ComplexRandomVariable and computes all the necessary 
 * values by numerical integration.  Might not be numerically accurate and is fairly slow.
 * THIS CURRENTLY ASSUMES THAT the true amplitude rho0 = 1
 */
class GeneralCLTComputer(override val M : Int, override val p : Double, val X : ComplexRandomVariable) extends AbstractCLTComputer(M,p) {
  
  //get the class representing the marginal magnitude random variable of X
  protected val Z = X.magnitudeMarginal
  
  override def fZ(z : Double) = Z.pdf(z)
  
  //Compute g by numerical integration
  override def g(phi : Double) = trapezoidal( r => r*f(r,phi), 0, 30*sqrt(Z.getVariance), 1000) 
  
  //Compute g2 by numerical integration
  override def g2(phi : Double) = trapezoidal( r => r*r*f(r,phi), 0, 30*sqrt(Z.getVariance), 1000) 
  
  /** The marginal pdf of the phase */
  def f(phi : Double) : Double = trapezoidal( r => f(r,phi), 0, 30*sqrt(Z.getVariance), 1000) 
  
}

/** 
 *CLT for circularly symmetric complex Gaussian noise.  Makes use of specific formula for
 *this case.  This is fast and should be numerically accurate.
*/
class GaussianCLT(override val M : Int, override val p : Double, val sigma : Double, val rho0 : Double) extends AbstractCLTComputer(M,p) {
  
  /** Signal amplitude defaults to 1 */
  def this(M : Int, p: Double, sigma : Double) = this(M,p,sigma,1.0)
  
  protected val k = rho0/sigma
  protected val Z = new Chi.Chi2(1/k/k)   //marginal distribution of the magnitude
  
  override def fZ(z : Double) = Z.pdf(z)
  
  //Compute g using known formula
  override def g(phi : Double) : Double = {
    //val a = cos(phi)
    //val mult = exp(-k*k/2)/(2*k)
    //val s1 = 2*a*k
    //val s2 = exp(sqr(a*k)/2)*sqrt(2*pi)*(1+sqr(a*k))*(1 + erf(a*k/sqrt(2)))
    //return mult*(s1+s2)/2/pi
    
    val a = cos(phi)  //equivalent formula from Barry
    val b = sin(phi)
    val PHI = (1 + erf(a*k/sqrt(2)))/2
    return a*exp(-k*k/2)/2/pi + 1/k/sqrt(2*pi) * exp(-k*k/2*b*b) * PHI * (1 + k*k*a*a)
  }
  
  //Compute g2 using known formula
  override def g2(phi : Double) : Double = {
    val a = cos(phi)
    val mult = exp(-k*k/2)/(2*k*k)
    val s1 = 4 + 2*sqr(a*k)
    val s2 = a*k*exp(sqr(a*k)/2)*sqrt(2*pi)*(3+sqr(a*k))*(1 + erf(a*k/sqrt(2)))
    return mult*(s1+s2)/2/pi
  }
  
  /** The marginal pdf of the phase */
  def f(phi : Double) : Double = {
    val a = cos(phi)
    val mult = exp(-k*k/2)/(2)
    val s1 = 2
    val s2 = a*k*exp(sqr(a*k)/2)*sqrt(2*pi)*(1 + erf(a*k/sqrt(2)))
    return mult*(s1+s2)/2/pi
  }
  
  //boost the accuracy of these integrals
  override lazy val A1 = trapezoidal( x => sqr(sin(x))*g2(x), -pi, pi, 50000)
  override lazy val A2 = trapezoidal( x => sqr(sin(fracpart(x)))*g2(x), -pi, pi, 50000)
  override def h2(x : Double) = trapezoidal( phi => cos(fracpart(x + phi))*g(phi), -pi, pi, 50000)
  override def h1(x : Double) = trapezoidal( phi => cos(x + phi)*g(phi), -pi, pi, 50000)
  override def m2(x : Double) = trapezoidal( phi => sin(fracpart(x + phi))*g(phi), -pi, pi, 50000)
  override def m1(x : Double) = trapezoidal( phi => sin(x + phi)*g(phi), -pi, pi, 50000)
  

}