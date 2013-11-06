disp('Building time offset estimator')
mex CXXFLAGS="\$CXXFLAGS -O3" LDXXFLAGS="\$LDXXFLAGS -O3" src/mexTimingEstimator.cpp src/MatlabHelper.cpp -Iinclude -outdir bin
%mex src/mexTimingEstimator.cpp src/MatlabHelper.cpp -Iinclude -outdir bin

disp('Compiling time offset estimator tests')
cd test/timeoffsettests
mex BrentTest.cpp -I../../include
mex FinitePulseTest.cpp -I../../include
mex RecursiveTest.cpp -I../../include
mex DirectTest.cpp -I../../include
mex UtilTest.cpp -I../../include
mex NaiveTest.cpp -I../../include
mex OerderMeyerMasseyTest.cpp -I../../include
cd oldtimingestimator
make_old_timingestimator
cd ..

disp('Running time offset estimator tests')
BrentTest
FinitePulseTest
RecursiveTest
DirectTest
UtilTest
NaiveTest
OerderMeyerMasseyTest
cd ..
mexTimingEstimator_test
cd ..