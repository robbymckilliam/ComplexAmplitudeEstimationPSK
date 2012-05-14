package cam.psk

import numbers.finite.Complex

/**
 * @author Robby McKilliam
 */
trait ComplexAmplitudeEstimator {
  
 def estimate(y : Seq[Complex] ) : Complex
 
 /** 
  * Returns the square error in amplitude and phase.
  */ 
 def error(ahat :Complex, a0 : Complex) : (Double, Double)
  
}

object ComplexAmplitudeEstimator {
 final def fracpart(x : Double) = x - scala.math.round(x)
}

/** 
 * Abstract class for Coherent PSK estimators.  Provides function for correctly computing
 * square error between phase and amplitude
 */
abstract class CoherentComplexAmplitudeEstimator extends ComplexAmplitudeEstimator {
  
  /** 
  * Returns the square error in amplitude and phase between two complex numbers.
  * Phase error is computed modulo 2pi.
  */ 
 final def error(ahat :Complex, a0 : Complex) : (Double, Double) = {
   val tpi = 2*scala.math.Pi
   val pe = tpi*ComplexAmplitudeEstimator.fracpart((ahat.angle - a0.angle)/tpi)
   val ae = ahat.magnitude - a0.magnitude
   return (ae*ae, pe*pe)
 }

  
}

/** 
 * Abstract class for Coherent PSK estimators.  Provides function for correctly computing
 * square error between phase and amplitude
 */
abstract class NonCoherentComplexAmplitudeEstimator(M : Int) extends ComplexAmplitudeEstimator {
  
  /** 
  * Returns the square error in amplitude and phase between two complex numbers.
  * Phase error is computed modulo 2pi/M.
  */ 
 final def error(ahat :Complex, a0 : Complex) : (Double, Double) = {
   val tpi = 2*scala.math.Pi/M
   val pe = tpi*ComplexAmplitudeEstimator.fracpart((ahat.angle - a0.angle)/tpi)
   val ae = ahat.magnitude - a0.magnitude
   return (ae*ae, pe*pe)
 }

}
