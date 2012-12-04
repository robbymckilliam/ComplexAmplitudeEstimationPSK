/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package cam.psk.uep

import org.junit.Test
import org.junit.Assert._

import numbers.finite.Complex
import numbers.finite.PolarComplex
import numbers.finite.RectComplex
import scala.math.Pi

class MackenthunUEPTest {

  val rand = new scala.util.Random
  val tol = 1e-6  
  
  @Test
  def testSameAsCoherentMackenthun() {   
    val a0 = new PolarComplex(1,rand.nextDouble*Pi*2)
    val L = 100
    val M = 4 //QPSK
    val s = (1 to L).map(m => new PolarComplex(1, 2*Pi*rand.nextInt(M)/M))
    
    //pilot position setup
    val numpilots = 10
    val P = 0 until numpilots //pilots at the front
    val D4 = numpilots until L //data at the back
    val D = Array(D4)
    val G = Array(M)
    val tester = new cam.psk.MackenthunCoherent(M,P,D4,s)
    val est = new MackenthunCoherentUEP(P,D,G,s)
    val iters = 100
    for(i <- 1 to iters) {
      val std = 0.01
      val y = s.map(si => a0*si + new RectComplex(std*rand.nextGaussian, std*rand.nextGaussian) ) 
      val ahat = est.estimate(y)
      val ahattester = tester.estimate(y)
      val err = (ahat - a0).mag2
      assertTrue((ahat-ahattester).magnitude < tol)
    }
  }
  
  @Test
  def testCoherentEstimate() {   
     val a0 = new PolarComplex(1,rand.nextDouble*Pi*2)
    val L = 100
    
    //pilot position setup
    val numpilots = 10
    val P = 0 until numpilots //pilots at the front
    val D2 = numpilots until L/2 //BPSK in the middle
    val D4 = L/2 until L //QPSK at the back
    val D = Array(D2, D4)
    val G = Array(2,4)
    
    val s = (0 until L).map{ i =>
      if(P.contains(i)) new PolarComplex(1, 2*Pi*rand.nextDouble) //random pilots
      else if( D2.contains(i) ) new PolarComplex(1, 2*Pi*rand.nextInt(2)/2) //generate a BPSK symbol
      else new PolarComplex(1, 2*Pi*rand.nextInt(4)/4) //generate a QPSK symbol
    }
    
    val est = new MackenthunCoherentUEP(P,D,G,s)
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
  
}
