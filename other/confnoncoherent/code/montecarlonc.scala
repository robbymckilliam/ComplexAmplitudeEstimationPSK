/**
* Run montecarlo.
*/
import cam.psk.MackenthunNonCoherent
import cam.psk.MackenthunCoherent
import cam.psk.ViterbiViterbi
import cam.psk.ComplexAmplitudeEstimator
import numbers.finite.RectComplex
import numbers.finite.Complex
import numbers.finite.PolarComplex
import cam.noise.ComplexGaussian
import pubsim.Util

val Ms = List(2,4) //BPSK, QPSK
val Ls = List(16, 256, 4096)
val a0 = new PolarComplex(1,2*scala.math.Pi*(new scala.util.Random).nextDouble)
val iters = 5000

//construct an array of noise distributions with a logarithmic scale
val SNRdBs = -15 to 20 by 1
val SNRs = SNRdBs.map(db => scala.math.pow(10.0, db/10.0))
val noises = SNRs.map( snr => new ComplexGaussian(a0.mag2/snr/2) ) //variance for real and imaginay parts (divide by 2)

val starttime = (new java.util.Date).getTime

println("Running noncoherent estimators simulation")
for( L <- Ls; M <- Ms ) {
    
    val D = 0 until L //data symbol indices

    //factory me
    def estfactory = List(
      () => new MackenthunNonCoherent(M,D)
      //() => new ViterbiViterbi( M, D, (d : Double) => 1.0 )
    )

    for( estf <- estfactory ) {
      
      val estname =  estf().getClass.getSimpleName + "M" + M + "L" + L.toString
      print("Running " + estname)
      val eststarttime = (new java.util.Date).getTime

      //for all the noise distributions (in parallel threads)
      val mselist = noises.par.map { noise =>	
	//generate some PSK symbols 
	val rand = new scala.util.Random
	val s = (1 to L).map(m => new PolarComplex(1, 2*scala.math.Pi*rand.nextInt(M)/M)) 
	val est = estf() //construct an estimator	  

	var msec = 0.0; var msea = 0.0; var msep = 0.0;
	for( itr <- 1 to iters ) {
	  //generate a recieved signal
	  val y = s.map(si => a0*si + noise.noise )
	  val ahat = est.estimate(y) //run the estimator 
	  val (ae, pe) = est.error(ahat, a0) //compute the error
	  msep += pe
	  msea += ae
	  msec += (ahat - a0).mag2 //compute the mean square error 
	}
	print(".")		      
	(msea/iters, msep/iters, msec/iters) //last thing is what gets returned
      }.toList
      
      val estruntime = (new java.util.Date).getTime - eststarttime
      println(" finished in " + (estruntime/1000.0) + " seconds.")
      
      val filea = new java.io.FileWriter("data/" + estname + "a")
      val filep = new java.io.FileWriter("data/" + estname + "p")
      val filec = new java.io.FileWriter("data/" + estname + "c")
      (mselist, SNRdBs).zipped.foreach{ (mse, snr) =>
	val (ma,mp,mc) = mse
				       filea.write(snr.toString.replace('E', 'e') + "\t" + ma.toString.replace('E', 'e')  + "\n") 
				       filep.write(snr.toString.replace('E', 'e') + "\t" + mp.toString.replace('E', 'e')  + "\n") 
				       filec.write(snr.toString.replace('E', 'e') + "\t" + mc.toString.replace('E', 'e')  + "\n") 
				     }
      filea.close; filep.close; filec.close; //close all the files we wrote to 

    }
}

//// Now run the unmodulated carrier
println
println("Running unmodulated carrier simulation")
for( L <- Ls; M <- Ms ) {
    
    val P = 0 until L //pilots at the front
    val D = L until L //data at the back

    //factory method to enable parallelism
    def estfactory = List(
      (p : Seq[Complex]) => new MackenthunCoherent(M,P,D,p)
    )

    for( estf <- estfactory ) {
      
      val estname =  estf(null).getClass.getSimpleName + "M" + M + "L" + L.toString + "absP" + P.length
      print("Running " + estname)
      val eststarttime = (new java.util.Date).getTime

      //for all the noise distributions (in parallel threads)
      val mselist = noises.par.map { noise =>	
	//generate some PSK symbols 
	val rand = new scala.util.Random
	val s = (1 to L).map(m => new PolarComplex(1, 2*scala.math.Pi*rand.nextInt(M)/M)) 
	val est = estf(s) //construct an estimator	  

	var msec = 0.0; var msea = 0.0; var msep = 0.0;
	for( itr <- 1 to iters ) {
	  //generate a recieved signal
	  val y = s.map(si => a0*si + noise.noise )
	  val ahat = est.estimate(y) //run the estimator 
	  val (ae, pe) = est.error(ahat, a0) //compute the error
	  msep += pe
	  msea += ae
	  msec += (ahat - a0).mag2 //compute the mean square error 
	}
	print(".")		      
	(msea/iters, msep/iters, msec/iters) //last thing is what gets returned
      }.toList
      
      val estruntime = (new java.util.Date).getTime - eststarttime
      println(" finished in " + (estruntime/1000.0) + " seconds.")
      
      val filea = new java.io.FileWriter("data/" + estname + "a")
      val filep = new java.io.FileWriter("data/" + estname + "p")
      val filec = new java.io.FileWriter("data/" + estname + "c")
      (mselist, SNRdBs).zipped.foreach{ (mse, snr) =>
	val (ma,mp,mc) = mse
				       filea.write(snr.toString.replace('E', 'e') + "\t" + ma.toString.replace('E', 'e')  + "\n") 
				       filep.write(snr.toString.replace('E', 'e') + "\t" + mp.toString.replace('E', 'e')  + "\n") 
				       filec.write(snr.toString.replace('E', 'e') + "\t" + mc.toString.replace('E', 'e')  + "\n") 
				     }
      filea.close; filep.close; filec.close; //close all the files we wrote to 

    }

}


val runtime = (new java.util.Date).getTime - starttime

println("Finished in " + (runtime/1000.0) + " seconds.\n") 





