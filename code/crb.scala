/**
  * Run montecarlo.
  */
import cam.psk.MackenthunCoherent
import cam.psk.ComplexAmplitudeEstimator
import cam.psk.BPSKCRB
import cam.psk.QPSKCRB
import cam.psk.MPSKCRB
import cam.psk.CramerRaoBound
import cam.noise.ComplexGaussian
import numbers.finite.RectComplex
import numbers.finite.Complex
import numbers.finite.PolarComplex
import pubsim.distributions.complex.ComplexRandomVariable
import pubsim.Util

/** Compute given CRB for given L and snrs */
def computecrb( crbf : (Double)=>CramerRaoBound, L : Int, snrdbs : Seq[Double], name : String) { 
  
  val SNRs = snrdbs.map(db => scala.math.pow(10.0, db/10.0))
  val vars = SNRs.map( snr => 1.0/snr/2 ) //variance for real and imaginay parts (divide by 2)

  print("Running " + name + " ")
  
  val crblist = vars.par.map { v =>
    val crbvar = crbf(v).variance(L)
			      print(".")
			      crbvar //last thing is what gets returned
			    }.toList

  val file = new java.io.FileWriter("data/" + name)
  (crblist, snrdbs).zipped.foreach{ (crb, snr) => file.write(snr.toString.replace('E', 'e') + "\t" + crb.toString.replace('E', 'e')  + "\n") }
  file.close
  println(" finished")

}

val Ls = List(32,256,2048,4096) //number of symbols
val SNRdBs = -20.0 to 20.0 by 0.5 //snrs we will use

//list of CRBS we will compute with names
val crbfs = List[((Double)=>CramerRaoBound, String)](
( (v : Double) => new BPSKCRB(v), "BPSKCRB" ),
( (v : Double) => new QPSKCRB(v), "QPSKCRB" ),
)

val starttime = (new java.util.Date).getTime

//run everything
for( L <- Ls; crbf <- crbfs ) computecrb(crbf._1, L, SNRdBs, crbf._2 + "" + L.toString)

//run a 8PSK with different SNR since it's a bit numerically sensitive at high SNR
computecrb( (v : Double) => new MPSKCRB(8,v), 4096, -5.0 to 15.0 by 0.5, "EightPSKCRB4096" )


val runtime = (new java.util.Date).getTime - starttime
println("Finished in " + (runtime/1000.0) + " seconds.\n")





