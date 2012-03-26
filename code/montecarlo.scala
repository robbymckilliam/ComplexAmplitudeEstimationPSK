/**
* Run montecarlo.
*/
import cam.psk.Mackenthun
import cam.psk.ComplexAmplitudeEstimator
import numbers.finite.RectComplex
import numbers.finite.Complex
import numbers.finite.PolarComplex
import pubsim.distributions.GaussianNoise
import pubsim.Util

val M = 4 //QPSK
val Ls = List(10, 100, 1000)
val a0 = new PolarComplex(1,2*scala.math.Pi*(new scala.util.Random).nextDouble)
val iters = 10000

//construct an array of noise distributions with a logarithmic scale
val SNRdBs = -10 to 20 by 1
val SNRs = SNRdBs.map(db => scala.math.pow(10.0, db/10.0))
val noises = SNRs.map( snr => new GaussianNoise(0.0, a0.mag2/snr/2) ) //variance for real and imaginay parts (divide by 2)

val starttime = (new java.util.Date).getTime
for( L <- Ls ) {

  //pilot position setup
  val numpilots = 5
  val P = 0 until numpilots //pilots at the front
  val D = numpilots until L //data at the back

  //factory me
  def estfactory = List(
    (p : Seq[Complex]) => new Mackenthun(M,P,D,p)
  )

  for( estf <- estfactory ) {
    
    val estname = estf(null).getClass.getSimpleName + "L" + L.toString + "absP" + numpilots
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
	val y = s.map(si => a0*si + new RectComplex(noise.getNoise, noise.getNoise) )
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
    
    val filea = new java.io.FileWriter(estname + "a")
    val filep = new java.io.FileWriter(estname + "p")
    val filec = new java.io.FileWriter(estname + "c")
    (mselist, SNRdBs).zipped.foreach{ (mse, snr) =>
	val (ma,mp,mc) = mse
	filea.write(snr.toString.replace('E', 'e') + "\t" + ma.toString.replace('E', 'e')  + "\n") 
	filep.write(snr.toString.replace('E', 'e') + "\t" + mp.toString.replace('E', 'e')  + "\n") 
        filec.write(snr.toString.replace('E', 'e') + "\t" + mc.toString.replace('E', 'e')  + "\n") 
    }
    filea.close; filep.close; filec.close //close all the files we wrote to 

    //println
    //for( i <- SNRdBs.indices ) println(SNRdBs(i) + "\t" + mselist(i))
  }

}

val runtime = (new java.util.Date).getTime - starttime

println("Finished in " + (runtime/1000.0) + " seconds.\n") 





