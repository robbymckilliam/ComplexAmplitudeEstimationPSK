% compile with internal messages as double (default)
mex CXXFLAGS="\$CXXFLAGS -O2" LDXXFLAGS="\$LDXXFLAGS -O2" ...
    mexOldTimeEstimator.cpp TimingEstimator.cpp MatlabHelper.cpp
% mex src/mexTimingEstimator.cpp src/TimingEstimator.cpp src/MatlabHelper.cpp -Iinclude -outdir bin
% compile with internal messages as float
% mex src/mexTimingEstimator.cpp src/TimingEstimator.cpp src/MatlabHelper.cpp -Iinclude -outdir bin -DVALTYPE=float -DMSGTYPE=float
