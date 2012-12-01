/*
 * @author Robby McKilliam
 */

package cam.noise

import org.junit._
import Assert._
import pubsim.distributions.complex.SymmetricComplexNormal

class CLTComputerTest {

  @Before
  def setUp: Unit = {
  }

  @After
  def tearDown: Unit = {
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
    for( phi <- 0.1 to 2*scala.math.Pi by 0.2) assertEquals(noise1.fZ(phi), noise2.fZ(phi), 0.0001)
    for( phi <- 0.1 to 2*scala.math.Pi by 0.2) assertEquals(noise1.g(phi), noise2.g(phi), 0.0001)
    for( phi <- 0.1 to 2*scala.math.Pi by 0.1) assertEquals(noise1.g2(phi), noise2.g2(phi), 0.0001)
    for( phi <- 0.1 to 2*scala.math.Pi by 0.2) assertEquals(noise1.f(phi), noise2.f(phi), 0.0001)
    for( phi <- 0.1 to 2*scala.math.Pi by 0.2) assertEquals(noise1.h1(phi), noise2.h1(phi), 0.0001)
    for( phi <- 0.1 to 2*scala.math.Pi by 0.2) assertEquals(noise1.h2(phi), noise2.h2(phi), 0.0001)
    for( phi <- 0.1 to 2*scala.math.Pi by 0.2) assertEquals(noise1.G(phi), noise2.G(phi), 0.0001)
    assertEquals(noise1.H, noise2.H, 0.0001)
    assertEquals(noise1.A1, noise2.A1, 0.0001)
    assertEquals(noise1.A2, noise2.A2, 0.0001)
    assertEquals(noise1.B1, noise2.B1, 0.0001)
    assertEquals(noise1.B2, noise2.B2, 0.0001)
    assertEquals(noise1.variance._1, noise2.variance._1, 0.001)
    assertEquals(noise1.variance._2, noise2.variance._2, 0.001)
   
    println(noise1.variance)
    println(noise2.variance)
  }
  
}