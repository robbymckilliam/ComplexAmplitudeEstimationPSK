#export JAVA_OPTS="-Xprof"
scala -nocompdaemon -cp lib/ComplexAmplitudePSK.jar:lib/ScalaNumbers.jar:lib/PubSim.jar:lib/Jama-1.0.2.jar:lib/flanagan.jar:lib/colt.jar:lib/RngPack.jar crb.scala
#scala -nocompdaemon -cp lib/ComplexAmplitudePSK.jar:lib/ScalaNumbers.jar:lib/PubSim.jar:lib/Jama-1.0.2.jar:lib/flanagan.jar:lib/colt.jar:lib/RngPack.jar montecarlonc.scala
