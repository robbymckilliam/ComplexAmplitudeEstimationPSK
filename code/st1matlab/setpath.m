% used to enable acces to sub-directory files
%--------------------------------------------------------------------------
% setpath.m - This script sets the paths for various modem modules.
% 
% setpath
%
% Inputs:
%       None
% Outputs:
%       None
% Comments:
%       This  script is called prior to each set of new bursts
%--------------------------------------------------------------------------
% Author:  ITR
% Project: ASRP
%--------------------------------------------------------------------------
% Copyright 2011 : 
% Institute for Telecommunications Research
% University of South Australia
%--------------------------------------------------------------------------
% Contains ITR Background IP
%--------------------------------------------------------------------------

rDIR = '.';

% addpath([rDIR '/SNR']);
% addpath([rDIR '/branches/st1demo_Flight/postproc']);
% addpath([rDIR '/branches/st1demo_Flight/postprocNTS']);
% addpath([rDIR '/branches/st1demo_Flight/visTools']);
addpath(rDIR);
addpath([rDIR '/.']);
addpath([rDIR '/Acquisition']);
addpath([rDIR '/impl/software/library/phy/bin']);
addpath([rDIR '/impl/software/ref-model/mlibrary']);
addpath([rDIR '/impl/software/ref-model/mlibrary/Acquisition']);
addpath([rDIR '/impl/software/ref-model/mlibrary/Acquisition/frequencyestimator']);
addpath([rDIR '/impl/software/ref-model/mlibrary/Acquisition/matchedfilter']);
addpath([rDIR '/impl/software/ref-model/mlibrary/Acquisition/noiseestimator']);
addpath([rDIR '/impl/software/ref-model/mlibrary/Acquisition/dcblockingfilter']);
addpath([rDIR '/impl/software/ref-model/mlibrary/Codecs/']);
clear rDIR
