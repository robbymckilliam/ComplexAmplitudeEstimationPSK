/**
* Run montecarlo.
*/
import cam.psk.uep.ViterbiViterbiUEP
import cam.psk.UnmodulatedCarrier
import cam.psk.MackenthunNonCoherent
import cam.psk.MackenthunCoherent
import cam.psk.uep.MackenthunCoherentUEP
import cam.psk.ViterbiViterbi
import cam.psk.ComplexAmplitudeEstimator
import numbers.finite.RectComplex
import numbers.finite.Complex
import numbers.finite.PolarComplex
import cam.noise.ComplexGaussian
import pubsim.Util
import scala.math.Pi

val Ls = List(30, 300, 3000)
val a0 = new PolarComplex(1,2*Pi*(new scala.util.Random).nextDouble)
val iters = 10000

//construct an array of noise distributions with a logarithmic scale
val SNRdBs = -20 to 20 by 1
val SNRs = SNRdBs.map(db => scala.math.pow(10.0, db/10.0))
val noises = SNRs.map( snr => new ComplexGaussian(a0.mag2/snr/2) ) //variance for real and imaginay parts (divide by 2)

val starttime = (new java.util.Date).getTime

println("Running noncoherent estimators simulation")
for( L <- Ls ) {
    
    val P = 0 until L/3 //pilot symbol indices
    val D2 = L/3 until 2*L/3 //BPSK symbol indices
    val D4 = 2*L/3 until L //BPSK symbol indices
    val D = Array(D2,D4)
    val G = Array(2,4)

    //factory method to enable parallelism
    def estfactory = List(
      (p : Seq[Complex]) => new UnmodulatedCarrier(0 until L,p),
      (p : Seq[Complex]) => new MackenthunCoherentUEP(P,D,G,p),
      (p : Seq[Complex]) => new ViterbiViterbiUEP(P,D,G,p, t => 1.0)
    )

    for( estf <- estfactory ) {
      
      val estname =  estf(null).getClass.getSimpleName + L.toString
      print("Running " + estname)
      val eststarttime = (new java.util.Date).getTime

      //for all the noise distributions (in parallel threads)
      val mselist = noises.par.map { noise =>	
	//generate some PSK symbols 
	val rand = new scala.util.Random
	val s = (0 until L).map{ i => 
	  if(P.contains(i)) new PolarComplex(1, 2*Pi*rand.nextDouble) //random pilots
	  else if( D2.contains(i) ) new PolarComplex(1, 2*Pi*rand.nextInt(2)/2) //generate a BPSK symbol
	  else new PolarComplex(1, 2*Pi*rand.nextInt(4)/4) //generate a QPSK symbol
	}
	val est = estf(s) //construct an estimator	  

	var msec = 0.0; var msea = 0.0; var msep = 0.0; var mseaunb = 0.0;
	for( itr <- 1 to iters ) {
	  //generate a recieved signal
	  val y = s.map(si => a0*si + noise.noise )
	  val ahat = est.estimate(y) //run the estimator 
	  val (ae, pe) = est.error(ahat, a0) //compute the error
	  msep += pe
	  msea += ae
	}
	print(".")		      
	(msea/iters, msep/iters) //last thing is what gets returned
      }.toList
      
      val estruntime = (new java.util.Date).getTime - eststarttime
      println(" finished in " + (estruntime/1000.0) + " seconds.")
      
      val filea = new java.io.FileWriter("data/" + estname + "a")
      val filep = new java.io.FileWriter("data/" + estname + "p")
      (mselist, SNRdBs).zipped.foreach{ (mse, snr) =>
	val (ma,mp) = mse
				       filea.write(snr.toString.replace('E', 'e') + "\t" + ma.toString.replace('E', 'e')  + "\n") 
				       filep.write(snr.toString.replace('E', 'e') + "\t" + mp.toString.replace('E', 'e')  + "\n") 
				     }
      filea.close; filep.close;//close all the files we wrote to 

    }
}

val runtime = (new java.util.Date).getTime - starttime
println("Finished in " + (runtime/1000.0) + " seconds.\n") 





