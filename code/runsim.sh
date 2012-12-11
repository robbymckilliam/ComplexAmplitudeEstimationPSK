#export JAVA_OPTS="-Xprof -server -Xms1g -Xmx1g"
export JAVA_OPTS="-d64 -server -Xms1g -Xmx1g"
scala -nocompdaemon -cp lib/ComplexAmplitudePSK.jar:lib/ScalaNumbers.jar:lib/PubSim.jar:lib/Jama-1.0.2.jar:lib/flanagan.jar:lib/colt.jar:lib/RngPack.jar:lib/jblas-1.2.0.jar:lib/junit-4.11.jar clt.scala
#scala -nocompdaemon -cp lib/ComplexAmplitudePSK.jar:lib/ScalaNumbers.jar:lib/PubSim.jar:lib/Jama-1.0.2.jar:lib/flanagan.jar:lib/colt.jar:lib/RngPack.jar:lib/jblas-1.2.0.jar:lib/junit-4.11.jar montecarlo.scala
#scala -nocompdaemon -cp lib/ComplexAmplitudePSK.jar:lib/ScalaNumbers.jar:lib/PubSim.jar:lib/Jama-1.0.2.jar:lib/flanagan.jar:lib/colt.jar:lib/RngPack.jar:lib/jblas-1.2.0.jar:lib/junit-4.11.jar montecarlonc.scala
