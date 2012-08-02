#export JAVA_OPTS="-Xprof"
scala -nocompdaemon -cp ../../../code/lib/ComplexAmplitudePSK.jar:../../../code/lib/ScalaNumbers.jar:../../../code/lib/PubSim.jar:../../../code/lib/Jama-1.0.2.jar:../../../code/lib/flanagan.jar:../../../code/lib/colt.jar:../../../code/lib/RngPack.jar montecarlonc.scala
