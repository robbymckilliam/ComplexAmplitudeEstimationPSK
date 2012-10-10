/*
 * @author Robby McKilliam
 */

package cam.noise

import org.junit._
import Assert._

class ComplexNoiseTest {

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
  def testComplexGaussian = {
    val M = 4
    val p = 0.1
    val rho0 = 1.0
    val sigma = 1.0
    val snr = rho0*rho0/2/sigma/sigma
    val kappa = rho0/sigma
    val noise = new ComplexGaussian(sigma*sigma)
    println(noise.magnitudeMarginal.getVariance)
    println(noise.noise)
  }
  
}
