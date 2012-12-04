/**
 * Mackenthun.scala
 * Implements various versions of Mackenthun's least squares carrier phase 
 * and amplitude estimator for phase-shift keyed signals.
 * @author Robby McKilliam
 */
package cam.psk

import numbers.finite.Complex
import numbers.finite.PolarComplex

/** 
 *  An implementation of Mackenthun's O(L log L) complex amplitude estimator.
 *  This is a minor modification of Mackenthun's estimator that allows for pilots.
 */
class MackenthunCoherent(val M: Int, val P: Seq[Int], val D : Seq[Int], val p : Seq[Complex]) 
extends CoherentComplexAmplitudeEstimator {
  
  protected val z = new Array[(Double, Int)](D.length)
  protected val g = new Array[Complex](D.length)
  
  protected val L = D.length + P.length
  
  protected val w = 2.0*scala.math.Pi/M
  protected val eta = new PolarComplex(1, w)
  protected val nu = eta - 1.0
  
  protected def sigma(k : Int)= z(k % D.length)._2

  /** Execute Mackenthun's estimator */
  def estimate(y : Seq[Complex]) : Complex = {
    
    //setup sequences and sort
    var Y = Complex.zero
    for(i <- P) Y = Y + y(i)*p(i).conjugate
    for( i <- D.indices ){
      val phi = y(D(i)).angle
      val u = w*scala.math.round(phi/w)
      z(i) = (phi - u,i)
      g(i) = y(D(i))*(new PolarComplex(1,-u))
      Y = Y + g(i)
    }
    scala.util.Sorting.quickSort(z) //scala sorts tuples in ascending order of the first element by default!
    
    //now run the search loop
    var ahat = Y / L
    var Qhat = Y.mag2 / L
    for( k <- 0 until M*D.length ){
      Y = Y + nu*g(sigma(k))
      g(sigma(k)) = g(sigma(k))*eta
      val Q = Y.mag2 / L
      if( Q > Qhat ){
        ahat = Y / L
        Qhat = Q
      }
    }
    
    return ahat
  }
  
}

/** 
 * Mackenthun's original NonCoherent estimator.
 * Does not allow for any pilot symbols.
 */
class MackenthunNonCoherent(override val M : Int, val D : Seq[Int]) 
extends NonCoherentComplexAmplitudeEstimator(M) {
  
  protected val z = new Array[(Double, Int)](D.length)
  protected val g = new Array[Complex](D.length)
  
  protected val L = D.length
  
  protected val w = 2.0*scala.math.Pi/M
  protected val eta = new PolarComplex(1, w)
  protected val nu = eta - 1.0
  
  protected def sigma(k : Int)= z(k % D.length)._2

  /** Execute Mackenthun's estimator */
  def estimate(y : Seq[Complex]) : Complex = {
    
    //setup sequences and sort
    var Y = Complex.zero
    for( i <- D.indices ){
      val phi = y(D(i)).angle
      val u = w*scala.math.round(phi/w)
      z(i) = (phi - u,i)
      g(i) = y(D(i))*(new PolarComplex(1,-u))
      Y = Y + g(i)
    }
    scala.util.Sorting.quickSort(z) //scala sorts tuples in ascending order of the first element by default!
    
    //now run the search loop
    var ahat = Y / L
    var Qhat = Y.mag2 / L
    for( k <- 0 until D.length ){
      Y = Y + nu*g(sigma(k))
      g(sigma(k)) = g(sigma(k))*eta
      val Q = Y.mag2 / L
      if( Q > Qhat ){
        ahat = Y / L
        Qhat = Q
      }
    }  
    
    return ahat
  }
}
