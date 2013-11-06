% compile with internal messages as double (default)
% mex src/mexQAMModulator.cpp src/QAMModulator.cpp src/MatlabHelper.cpp -Iinclude -outdir bin
% compile with internal messages as float
% mex src/mexQAMModulator.cpp src/QAMModulator.cpp src/MatlabHelper.cpp -Iinclude -outdir bin -DMSGTYPE=float

mex -output mexQAMModulator src/mexQAMModulatorV2.cpp src/QAMModulatorV2.cpp src/MatlabHelper.cpp -Iinclude -outdir bin -DMSGTYPE=float

cd bin
mexQAMModulator
cd ..

cd test
test_QAMModulator;
cd ..
