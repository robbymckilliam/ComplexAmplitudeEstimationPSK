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
  def testNoise= {
    val M = 4
    val p = 0.1
    val rho0 = 1.0
    val sigma = 1.0
    val snr = rho0*rho0/2/sigma/sigma
    val kappa = rho0/sigma
    val X = new SymmetricComplexNormal(sigma*sigma)
    val noise = new GeneralCLTComputer(M,p,X)
    println("H " + noise.H)
    println("A1 " + noise.A1)
    println("A2 " + noise.A2)
    println("g")
    for( x <- 0.0 to 2.0*scala.math.Pi by 0.2 ) println(noise.g(x))
    //println("G")
    //for( x <- 0.0 to 2.0*scala.math.Pi by 2.0 ) println(noise.G(x))
  }
  
}