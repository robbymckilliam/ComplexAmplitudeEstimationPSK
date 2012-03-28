
package cam.psk

import org.junit.Test
import org.junit.Assert._

import numbers.finite.Complex
import numbers.finite.PolarComplex
import numbers.finite.RectComplex

/** 
 * @author Robby McKilliam
 */
class MackenthunTest {

  val rand = new scala.util.Random
  
  val a0 = new PolarComplex(1,rand.nextDouble*scala.math.Pi*2)
  val L = 100
  val M = 4 //QPSK
  val s = (1 to L).map(m => new PolarComplex(1, 2*scala.math.Pi*rand.nextInt(M)/M))
    
  //pilot position setup
  val numpilots = 10
  val P = 0 until numpilots //pilots at the front
  val D = numpilots until L //data at the back
  
  @Test
  def testCoherentEstimate() {   
    val est = new MackenthunCoherent(M,P,D,s)
    val iters = 100
    var mse = 0.0
    for(i <- 1 to iters) {
      val std = 0.01
      val y = s.map(si => a0*si + new RectComplex(std*rand.nextGaussian, std*rand.nextGaussian) ) 
      val ahat = est.estimate(y)
      //println(ahat + " and " + a0)
      val err = (ahat - a0).mag2
      mse = mse + err
      assertTrue(err < 0.01)
    }
    println(mse/iters)
  }
  
  @Test
  def testNonCoherentEstimate() {      
    val estnc = new MackenthunNonCoherent(M,D)
    val iters = 100
    var mse = 0.0
    for(i <- 1 to iters) {
      val std = 0.01
      val y = s.map(si => a0*si + new RectComplex(std*rand.nextGaussian, std*rand.nextGaussian) ) 
      val ahat = estnc.estimate(y)
      //println(ahat + " and " + a0)
      val (ae, pe) = estnc.error(ahat, a0)
      mse = mse + pe
      assertTrue(pe < 0.01)
    }
    println(mse/iters)
  }
  
  
}
