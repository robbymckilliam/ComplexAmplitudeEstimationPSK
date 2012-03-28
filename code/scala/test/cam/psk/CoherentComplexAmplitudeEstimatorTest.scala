/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package cam.psk

import numbers.finite.PolarComplex

import org.junit._
import Assert._

class CoherentComplexAmplitudeEstimatorTest {

  @Before
  def setUp: Unit = {
  }

  @After
  def tearDown: Unit = {
  }

  @Test
  def example = {
  }
  
  @Test
  def testError = {
    val M = 4
    val P = List(); val D = List();
    val est = new MackenthunCoherent(M,P,D,null)
    
    {
      val ahat = new PolarComplex(1.0, 0.1)
      val a0 = new PolarComplex(1.1, -0.1)
      val (ae, pe) = est.error(ahat, a0)
      assertEquals(ae, 0.1*0.1, 0.00000001)
      assertEquals(pe, 0.2*0.2, 0.00000001)
    }
    
    {
      val ahat = new PolarComplex(1.0, scala.math.Pi + 0.1)
      val a0 = new PolarComplex(2.0, scala.math.Pi - 0.1)
      val (ae, pe) = est.error(ahat, a0)
      assertEquals(ae, 1.0, 0.00000001)
      assertEquals(pe, 0.2*0.2, 0.00000001)
    }
    
    {
      val ahat = new PolarComplex(1.0, -scala.math.Pi + 0.1)
      val a0 = new PolarComplex(2.0, scala.math.Pi - 0.1)
      val (ae, pe) = est.error(ahat, a0)
      assertEquals(ae, 1.0, 0.00000001)
      assertEquals(pe, 0.2*0.2, 0.00000001)
    }
    
  }

}
