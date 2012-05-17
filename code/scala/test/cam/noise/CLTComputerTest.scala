/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package cam.noise

import org.junit._
import Assert._
import pubsim.distributions.complex.ComplexRandomVariable
import pubsim.distributions.complex.SymmetricComplexNormal

class CLTComputerTest {

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
  def testGaussianCLT = {
    val M = 4
    val p = 0.1
    val rho0 = 1.0
    val sigma = 1.0
    val snr = rho0*rho0/2/sigma/sigma
    val kappa = rho0/sigma
    val X = new SymmetricComplexNormal(sigma*sigma)
    val noise1 = new GeneralCLTComputer(M,p,X)
    val noise2 = new GaussianCLT(M,p,sigma)
    assertEquals(noise1.fZ(0.1), noise2.fZ(0.1), 0.001)
    assertEquals(noise1.g(0.1), noise2.g(0.1), 0.001)
    assertEquals(noise1.g2(0.1), noise2.g2(0.1), 0.001)
    assertEquals(noise1.H, noise2.H, 0.001)
    assertEquals(noise1.A1, noise2.A1, 0.001)
    assertEquals(noise1.A2, noise2.A2, 0.001)
  }
  
}