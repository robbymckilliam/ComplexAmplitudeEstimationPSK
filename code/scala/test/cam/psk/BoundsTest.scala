/*
 * Test for some of the bounds
 */

package cam.psk

import org.junit._
import Assert._
import numbers.finite.RectComplex
import scala.math.sqrt

class BoundsTest {

  @Before
  def setUp: Unit = {
  }

  @After
  def tearDown: Unit = {
  }

  @Test
  def testFirstQuadrant = {
    val q1 = cam.psk.MPSKCRB.firstquadrant(4)
    assertEquals(q1.length,1)
    assertTrue( (q1(0) - RectComplex(sqrt(2)/2,sqrt(2)/2)).mag2 < 0.000001 )
  }

}
