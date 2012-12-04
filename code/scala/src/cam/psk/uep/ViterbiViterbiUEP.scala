/*
 * @author Robby McKilliam
 */

package cam.psk.uep

import scala.math.Pi
import scala.math.round
import numbers.finite.PolarComplex
import numbers.finite.Complex
import cam.psk.CoherentComplexAmplitudeEstimator

class ViterbiViterbiUEP(val P: Seq[Int], val D : Seq[Seq[Int]], val G : Seq[Int], val p : Seq[Complex], val F : Double => Double) 
extends CoherentComplexAmplitudeEstimator {
  
  val absP = P.length //number of pilot symbols
  val absD = D.foldLeft(0)((s,e) => s + e.length) //number of data symbols
  val L = absP + absD //the total number of symbols
  
  val K = if( P.size > 0 ) 1.0 else G.min //scale factor for this estimator (see document) 
  
  def estimate(y : Seq[Complex] ) : Complex = {
    var A = Complex.zero
    for( i <- P ) A = A + y(i)*p(i).conjugate/y(i).magnitude*F(y(i).magnitude)
    for( ell <- G.indices; i <- D(ell) ) {
      val m = G(ell)
      A = A + (y(i)/y(i).magnitude).pow(m) * F(y(i).magnitude)
    }
    return new PolarComplex(A.magnitude, A.angle / K)
  }

}
