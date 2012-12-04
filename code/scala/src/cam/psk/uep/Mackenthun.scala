/*
 * @author Robby McKilliam
 */

package cam.psk.uep

import scala.math.Pi
import scala.math.round
import numbers.finite.PolarComplex
import numbers.finite.Complex
import cam.psk.CoherentComplexAmplitudeEstimator

/**
 * A modification of Mackenthun's estimator to allow for pilot symbols and data symbols from
 * multiple phase shift keyed constellations.  See the document:
 * Robby McKilliam, Andre Pollok, Bill Cowley "Carrier phase and amplitude estimation for phase
 * shift keying with unequal error protection"
 * The notation here follows that in the document, expect that sequences are indexed beginning
 * at 0 rather than 1 and the definition of D and G are slightly different.  Here D(i) is a sequence of 
 * indicies modulated using a constellation of size G(i).
 */
class MackenthunCoherentUEP(val P: Seq[Int], val D : Seq[Seq[Int]], val G : Seq[Int], val p : Seq[Complex])
extends CoherentComplexAmplitudeEstimator {
  
  val absP = P.length //number of pilot symbols
  val absD = D.foldLeft(0)((s,e) => s + e.length) //number of data symbols
  val H = G.indices.foldLeft(0)( (s,i) => s + G(i)*D(i).length ) //the size of the sort (see H in the paper)
  val L = absP + absD //the total number of symbols

  protected val t = new Array[(Double,Int,Int)](H) //the array of triples
  protected val g = new Array[Complex](L)
  
  final def m(k : Int) = t(k)._3
  final def sigma(k : Int) = t(k)._2
  final def roundm(x : Double, m : Int) = 2*Pi/m * round(x*m/(2*Pi))
  
  def estimate(y : Seq[Complex] ) : Complex = {
    
    //setup sequences and sort
    var Y = Complex.zero
    var c = 0
    for(i <- P) Y = Y + y(i)*p(i).conjugate
    for( ell <- G.indices ) {
      for( i <- D(ell) ) {
        val m = G(ell)
        val phi = y(i).angle
        val u = roundm(phi,m)
        g(i) = y(i)*(new PolarComplex(1,-u))
        Y = Y + g(i)
        val z = phi - u
        for( k <- 0 until m ) {
          t(c) = (z + 2*Pi*(k+1)/m, i, m) //store the triples
          c = c + 1
        }
      }
    }
    scala.util.Sorting.quickSort(t) //scala sorts tuples in ascending order of the first element by default!
    
    //now run the search loop
    var ahat = Y / L
    var Qhat = Y.mag2 / L
    for( k <- 0 until H ){
      val eta = new PolarComplex(1,2*Pi/m(k))
      Y = Y + (eta-1)*g(sigma(k))
      g(sigma(k)) = eta*g(sigma(k))
      val Q = Y.mag2 / L
      if( Q > Qhat ){
        ahat = Y / L
        Qhat = Q
      }
    }
    
    return ahat
  }
  
}
