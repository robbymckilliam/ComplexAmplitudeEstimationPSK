#g++ -pg profile.cpp ./src/CoherentMackenthun.cpp ./src/IndexedReal.cpp -Iinclude -o profile
#chmod 777 profile
#./profile
#gprof profile gmon.out > analysis.txt

g++ -std=c++0x -pthread -O3 turbosyncsim.cpp -IC -o runturbo
chmod 777 runturbo
./runturbo