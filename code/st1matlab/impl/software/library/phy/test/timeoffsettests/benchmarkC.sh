g++ -pg -std=c++0x benchmark.cpp -I../../include -o benchmark
chmod 777 benchmark
./benchmark
gprof benchmark gmon.out > analysis.txt

#g++ -std=c++0x -O3 benchmark.cpp -I../../include -o benchmark
#chmod 777 benchmark
#./benchmark