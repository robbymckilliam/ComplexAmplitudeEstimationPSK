#export JAVA_OPTS="-Xprof -d64 -server -Xms1g -Xmx1g"
export JAVA_OPTS="-d64 -server -Xms1g -Xmx1g"
scala -nocompdaemon -cp ../../../code/lib/ComplexAmplitudePSK.jar:../../../code/lib/ScalaNumbers.jar:../../../code/lib/PubSim.jar:../../../code/lib/Jama-1.0.2.jar:../../../code/lib/flanagan.jar:../../../code/lib/colt.jar:../../../code/lib/RngPack.jar:../../../code/lib/jblas-1.2.0.jar:../../../code/lib/junit-4.11.jar montecarlo.scala
