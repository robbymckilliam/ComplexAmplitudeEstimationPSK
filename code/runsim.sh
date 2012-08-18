#export JAVA_OPTS="-Xprof"
scala -nocompdaemon -cp lib/ComplexAmplitudePSK.jar:lib/ScalaNumbers.jar:lib/PubSim.jar:lib/Jama-1.0.2.jar:lib/flanagan.jar:lib/colt.jar:lib/RngPack.jar clt.scala
scala -nocompdaemon -cp lib/ComplexAmplitudePSK.jar:lib/ScalaNumbers.jar:lib/PubSim.jar:lib/Jama-1.0.2.jar:lib/flanagan.jar:lib/colt.jar:lib/RngPack.jar montecarlo.scala
scala -nocompdaemon -cp lib/ComplexAmplitudePSK.jar:lib/ScalaNumbers.jar:lib/PubSim.jar:lib/Jama-1.0.2.jar:lib/flanagan.jar:lib/colt.jar:lib/RngPack.jar montecarlonc.scala
