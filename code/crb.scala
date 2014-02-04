/**
  * Run montecarlo.
  */
import cam.psk.MackenthunCoherent
import cam.psk.ComplexAmplitudeEstimator
import cam.psk.BPSKCRB
import cam.noise.ComplexGaussian
import numbers.finite.RectComplex
import numbers.finite.Complex
import numbers.finite.PolarComplex
import pubsim.distributions.complex.ComplexRandomVariable
import pubsim.Util

val Ls = List(32,256,2048,4096)

//construct an array of noise distributions with a logarithmic scale
val SNRdBs = -20.0 to 20.0 by 0.5
val SNRs = SNRdBs.map(db => scala.math.pow(10.0, db/10.0))
val vars = SNRs.map( snr => 1.0/snr/2 ) //variance for real and imaginay parts (divide by 2)

val starttime = (new java.util.Date).getTime

for( L <- Ls ) {

  { //BPSK first
    print("Computing BPSK CRBs L = " + L + " ")
    
    val crblist = vars.par.map { v =>
      val crbvar = (new BPSKCRB).variance(v,L)
      print(".")
      crbvar //last thing is what gets returned
    }.toList

    val file = new java.io.FileWriter("data/BPSKCRB" + L)
      (crblist, SNRdBs).zipped.foreach{ (crb, snr) =>
        file.write(snr.toString.replace('E', 'e') + "\t" + crb.toString.replace('E', 'e')  + "\n")
      }
    file.close;
    println(" finished")
  }

}


val runtime = (new java.util.Date).getTime - starttime

println("Finished in " + (runtime/1000.0) + " seconds.\n")





