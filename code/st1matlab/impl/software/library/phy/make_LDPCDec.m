% compile with internal messages as double (default)
mex src/mexLDPCDec.cpp src/LDPCDec.cpp src/MatlabHelper.cpp -Iinclude -outdir bin
% compile with internal messages as float
% mex src/mexLDPCDec.cpp src/LDPCDec.cpp src/MatlabHelper.cpp -Iinclude -outdir bin -DMSGTYPE=float

cd bin
mexLDPCDec
cd ..

cd test
test_LDPCDec;
cd ..
