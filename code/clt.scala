/**
* Run montecarlo.
*/
import cam.psk.MackenthunCoherent
import cam.psk.ComplexAmplitudeEstimator
import cam.noise.ComplexGaussian
import numbers.finite.RectComplex
import numbers.finite.Complex
import numbers.finite.PolarComplex
import pubsim.distributions.complex.ComplexRandomVariable
import pubsim.Util

val Ms = List(2,4,8) //BPSK, QPSK, 8-PSK
//val Ls = List(2,4,8,16,32,64,128,256,512,1024,2048)
//val Ls = List(32,128,256,1024,4096)
val Ls = List(128)
val a0 = new PolarComplex(1,2*scala.math.Pi*(new scala.util.Random).nextDouble)

//construct an array of noise distributions with a logarithmic scale
val SNRdBs = -20.0 to 20.0 by 0.5
val SNRs = SNRdBs.map(db => scala.math.pow(10.0, db/10.0))
val noises = SNRs.map( snr => new ComplexGaussian(a0.mag2/snr/2) ) //variance for real and imaginay parts (divide by 2)

val starttime = (new java.util.Date).getTime
for( L <- Ls; M <- Ms ) {

  //for a range of different numbers of pilots
  for( numpilots <- List( 0, L/64, L/32, L/16, L/8, L/4, L/2, L ).distinct ) {

    val P = 0 until numpilots //pilots at the front
    val D = numpilots until L //data at the back

    //now compute the clts
    val cltname = "clt" + "M" + M + "L" + L + "absP" + P.length
    print("Computing " + cltname)
    //for all the noise distributions (in parallel threads)
    val cltlist = noises.par.map { noise =>	
	val p = P.length.toDouble/L  
	val (varp, vara) : (Double, Double) = noise.clt(M,p).variance	      
	print(".")
	(varp/L, vara/L) //last thing is what gets returned
      }.toList

     val filep = new java.io.FileWriter("data/" + cltname + "p")
     val filea = new java.io.FileWriter("data/" + cltname + "a")
     (cltlist, SNRdBs).zipped.foreach{ (mse, snr) =>
	val (vp,va) = mse
				       filea.write(snr.toString.replace('E', 'e') + "\t" + va.toString.replace('E', 'e')  + "\n") 
				       filep.write(snr.toString.replace('E', 'e') + "\t" + vp.toString.replace('E', 'e')  + "\n") 
				     }
     filea.close; filep.close //close all the files we wrote to 
    println(" finished")

  }

}

val runtime = (new java.util.Date).getTime - starttime

println("Finished in " + (runtime/1000.0) + " seconds.\n") 





