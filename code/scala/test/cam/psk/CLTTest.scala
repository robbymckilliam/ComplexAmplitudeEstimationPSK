/*
 * Test for the CLT computers
 */

package cam.psk

import org.junit._
import Assert._
import numbers.finite.PolarComplex
import pubsim.distributions.complex.ComplexRandomVariable
import pubsim.distributions.complex.SymmetricComplexNormal


class CLTTest {

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
  def testCoherentCLT = {
    val M = 4
    val p = 0.1
    val rho0 = 1.0
    val sigma = 0.1
    val snr = rho0*rho0/2/sigma/sigma
    val kappa = rho0/sigma
    val X = new SymmetricComplexNormal(sigma*sigma)
    val clt = new CoherentCLT(M,X,p)
    println(clt.variance)
  }

}
