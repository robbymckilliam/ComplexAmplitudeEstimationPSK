#g++ -pg profile.cpp ./src/CoherentMackenthun.cpp ./src/IndexedReal.cpp -Iinclude -o profile
#chmod 777 profile
#./profile
#gprof profile gmon.out > analysis.txt

g++ -std=c++11 -O3 turbosyncsim.cpp C/LDPCDec.cpp C/CoherentMackenthun.cpp -pthread -Wl,--no-as-needed -IC -o runturbo
chmod 777 runturbo
./runturbo
