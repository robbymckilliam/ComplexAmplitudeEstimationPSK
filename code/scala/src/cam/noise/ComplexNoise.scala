  /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package cam.noise

import numbers.finite.Complex
import numbers.finite.RectComplex

trait ComplexNoise extends pubsim.distributions.complex.ComplexRandomVariable {
  
  //Convert the complex number that comes from X, a ComplexRandomVariable
  def noise : Complex = {
    val c : pubsim.Complex = this.getNoise
    return new RectComplex(c.re,c.im)
  }
  
  /** 
   * Return a class for computing function relavant for the asymptotic distritbution of the least
   * estimator.  Requires M (for M-PSK) and also p, the proportion of pilots.
  */
  def clt(M : Int, p : Double) : CLTComputer
  
}

class ComplexGaussian(variance : Double) 
  extends pubsim.distributions.complex.SymmetricComplexNormal(variance) with ComplexNoise {
    
  /** Assumes signal amplitude rho0 is 1 */
  override def clt(M : Int, p : Double) : CLTComputer = {
    val sigma = scala.math.sqrt(variance)
    return new GaussianCLT(M,p,sigma)
  }
  
}
