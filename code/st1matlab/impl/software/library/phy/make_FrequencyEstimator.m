%% compile with double precision (default)
mex CXXFLAGS="\$CXXFLAGS -Wall -O3 -I/usr/local/include" LDFLAGS="\$LDFLAGS /usr/local/lib/libfftw3.a -lm" ...
    src/mexFrequencyEstimator.cpp src/FrequencyEstimator.cpp src/MatlabHelper.cpp -Iinclude -outdir bin

%% compile with single precision (float)
% CAUTION: the mex interface is currently buggy when compiled with float:
%          the C code returns the correct complex gain estimate, but the
%          version returned to Matlab is incorrect!
% mex CXXFLAGS="\$CXXFLAGS -Wall -O3 -I/usr/local/include" LDFLAGS="\$LDFLAGS /usr/local/lib/libfftw3f.a -lm" ...
%     src/mexFrequencyEstimator.cpp src/FrequencyEstimator.cpp src/MatlabHelper.cpp -Iinclude -outdir bin -DVALTYPE=float -DMSGTYPE=float -DFFTWFLOAT=1

%% test estimator
cd bin
mexFrequencyEstimator
cd ..

cd test
addpath('../bin/');
FrequencyEstimatorTest;
cd ..
