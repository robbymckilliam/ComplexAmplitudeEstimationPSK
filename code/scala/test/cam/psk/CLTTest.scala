/*
 * Test for the CLT computers
 */

package cam.psk

import org.junit._
import Assert._
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
    val sigma = 1.0
    val snr = rho0*rho0/2/sigma/sigma
    val kappa = rho0/sigma
    val X = new SymmetricComplexNormal(sigma*sigma)
    val clt = new CoherentCLT(M,X,p)
    println("H " + clt.H)
    println("A1 " + clt.A1)
    println("A2 " + clt.A2)
    println("var " + clt.variance)
    println("Z var " + clt.Z.getVariance)
    println("g")
    for( x <- 0.0 to 2.0*scala.math.Pi by 0.2 ) println(clt.g(x))
    println()
    println("fZ")
    for( z <- 0.0 to 1000*clt.Z.getVariance by 10*clt.Z.getVariance ) println(clt.fZ(z))
    //println()
    //println("f")
    //for( z <- 0.0 to 10.0 by 0.2 ) println(clt.f(z,0.0))
    
  }

}
