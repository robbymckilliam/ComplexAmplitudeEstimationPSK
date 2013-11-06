% compile with internal messages as double (default)
%mex CXXFLAGS="\$CXXFLAGS -std=c++0x" LDXXFLAGS="\$LDXXFLAGS -std=c++0x" ...
%    src/mexCoherentMackenthun.cpp src/IndexedReal.cpp
%    src/CoherentMackenthun.cpp src/MatlabHelper.cpp -Iinclude
%    -outdir bin
mex CXXFLAGS="-Wall \$CXXFLAGS -O2" LDXXFLAGS="\$LDXXFLAGS -O2" ...
    src/mexCoherentMackenthun.cpp src/IndexedReal.cpp src/CoherentMackenthun.cpp src/MatlabHelper.cpp -Iinclude -outdir bin
%mex src/mexCoherentMackenthun.cpp src/IndexedReal.cpp src/CoherentMackenthun.cpp src/MatlabHelper.cpp -Iinclude -outdir bin
% compile with internal messages as float
% mex src/mexCoherentMackenthun.cpp srcCoherentMackenthun.cpp src/MatlabHelper.cpp -Iinclude -outdir bin -DMSGTYPE=float

cd bin
%run binary without parameters (just displays help)
mexCoherentMackenthun
cd ..

cd test
disp('Compiling tests')
%mex CXXFLAGS="\$CXXFLAGS -std=c++0x" LDXXFLAGS="\$LDXXFLAGS -std=c++0x" mextest.cpp ../src/CoherentMackenthun.cpp ../src/MatlabHelper.cpp ...
%    -I../include
mex mextest.cpp ../src/CoherentMackenthun.cpp ../src/MatlabHelper.cpp ../src/IndexedReal.cpp -I../include

%run mex tests
disp('Running C tests')
mextest

%run matlab tests
disp('Running matlab tests')
matlabtests

cd ..
