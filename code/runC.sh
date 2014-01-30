#g++ -pg profile.cpp ./src/CoherentMackenthun.cpp ./src/IndexedReal.cpp -Iinclude -o profile
#chmod 777 profile
#./profile
#gprof profile gmon.out > analysis.txt

g++ -std=c++0x -pthread -O3 montecarlo.cpp -IC -o montecarlo
chmod 777 montecarlo
./montecarlo