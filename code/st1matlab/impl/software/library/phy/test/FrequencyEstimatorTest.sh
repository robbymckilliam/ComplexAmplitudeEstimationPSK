##############################################
#   NO PROFILER
##############################################
# compile with double precision
g++ -O3 FrequencyEstimatorTest.cpp ../src/FrequencyEstimator.cpp -I../include -o ../bin/FrequencyEstimatorTest -lfftw3
# compile with single precision
#g++ -O3 FrequencyEstimatorTest.cpp ../src/FrequencyEstimator.cpp -I../include -o ../bin/FrequencyEstimatorTest -lfftw3f -DVALTYPE=float -DMSGTYPE=float -DFFTWFLOAT=1
chmod 777 ../bin/FrequencyEstimatorTest
../bin/FrequencyEstimatorTest

##############################################
#   WITH PROFILER
##############################################
# compile with double precision & run profiler
#g++ -O3 FrequencyEstimatorTest.cpp ../src/FrequencyEstimator.cpp -I../include -o ../bin/FrequencyEstimatorTestProf -lfftw3 -pg
# compile with single precision & run profiler
#g++ -O3 FrequencyEstimatorTest.cpp ../src/FrequencyEstimator.cpp -I../include -o ../bin/FrequencyEstimatorTestProf -lfftw3f -DVALTYPE=float -DMSGTYPE=float -DFFTWFLOAT=1 -pg
#chmod 777 ../bin/FrequencyEstimatorTestProf
#../bin/FrequencyEstimatorTestProf
#gprof --flat-profile ../bin/FrequencyEstimatorTestProf gmon.out > ./profFrequencyEstimator.txt