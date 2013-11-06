function ST1ULacquisition_init(acqparams)
% This function initialises the ST1UL acquisition module and must be called
% prior to running pskacquisition.m. Specifically, this function constructs
% objects for the time and phase offset estimators.
%
% Comment: It is possible to only use a subset of pilot and data symbols in
%          the timing and phase estimators. For example, this can be useful
%          to meet constraints on the computational complexity.
%
%   ST1ULacquisition_init(acqparams)
%
% INPUTS:
%   acqparams       - struct holding aquisition parameters
%       .T              - symbol period
%       .OS             - oversampling factor
%       .P              - position (indices) of pilot symbols
%       .p              - pilot symbols
%       .D              - position (indices) of data symbols
%       .PTime          - position (indices) of pilot symbols used by
%                         timing estimator (equal to or a subset of .P)
%       .pTime          - pilot symbols used by timing estimator (equal to 
%                         or a subset of .p)
%       .DTime          - position (indices) of data symbols used by timing
%                         estimator (equal to or a subset of .D)
%       .PPhase         - position (indices) of pilot symbols used by phase
%                         estimator (equal to or a subset of .P)
%       .pPhase         - pilot symbols used by phase estimator (equal to
%                         or a subset of .p)
%       .DPhase         - position (indices) of data symbols used by phase
%                         estimator (equal to or a subset of .D)
%       .M              - M for M-PSK, i.e. size of constellation
%       .rolloff        - roll-off factor of root raised cosine pulse
%       .Tpshape        - duration of the transmit pulse shape in seconds
%       .timingEstimatorName
%                       - string identifying the timing estimator
%       .taumin         - minimum time offset
%       .taumax         - maximum time offset
%       .phaseEstimatorName
%                       - string identifying the phase estimator
%       .fr_max         - maximum Doppler rate
%       .fr_min         - minimum Doppler rate
%       .f_max          - range of residual Doppler after initial frequency
%                         estimation
%       .f_min          - range of residual Doppler after initial frequency
%                         estimation
%       .searchOS       - oversampling factor determining the search grid
%                         point spacing in the phase estimator. Need not be
%                         integer. High vlaues give higher accuracy, but
%                         also higher computational complexity. A good
%                         starting point is 2.
%       .Nfft           - FFT size for initial frequency estimator
%
% OUTPUTS:
%
% (c) 2012 ITR-UniSA
% Authors: Andre Pollok
% Created: 1 August 2012
% Last Modified: 5 March 2013

%% Map input QPSK constellation to internal constellation
% This function internally uses a QPSK constellation where the
% constellation points lie on the axes. In contrast Marc's framework uses a
% QPSK constellation where the points are halfway between the axes, i.e.
% rotated by 45 degrees. The mismatch is compensated by rotating the pilot
% symbols by 45 degrees. Furthermore, before returning the matched filter
% output and complex gain estimate, they are adjusted to match the input
% constellation.

if (acqparams.M == 4)
    % rotate pilot symbols by 45 degrees
    acqparams.p      = acqparams.p      * exp(1i*pi/4);
    acqparams.pTime  = acqparams.pTime  * exp(1i*pi/4);
    acqparams.pPhase = acqparams.pPhase * exp(1i*pi/4);
end

%% construct timing estimator object

T = acqparams.T; % symbol period
OS = acqparams.OS; % oversampling factor
rolloff = acqparams.rolloff; % roll-off factor
Tpshape = acqparams.Tpshape; % pulse duration

timingEstimatorName = acqparams.timingEstimatorName;
if ~strcmpi(timingEstimatorName,'disable')
    % pilots and data symbol indices
    DTime = acqparams.DTime;
    PTime = acqparams.PTime;
    pTime = acqparams.pTime;
    % adjuct search window:
    %   The timing estimator assumes the first sample at t=Ts and the first
    %   symbol at t=T=OS*Ts, i.e. it does not take the leading Tx filter tail
    %   into account. Therefore, we need to shift the search window by the
    %   leading Tx filter tail minus (OS-1) sample periods.
    taucorrection = acqparams.Tpshape/2 - (OS-1)*T/OS;
    taumin = acqparams.taumin + taucorrection;
    taumax = acqparams.taumax + taucorrection;
    
    mexTimingEstimator(timingEstimatorName,PTime,DTime,pTime,T,T/OS,1,OS,taumin,taumax,rolloff,Tpshape);
end

%% construct phase estimator object

phaseEstimatorName = acqparams.phaseEstimatorName;
if ~strcmpi(phaseEstimatorName,'disable') && ~strcmpi(phaseEstimatorName,'data-aided')
    % pilots and data symbol indices
    DPhase = acqparams.DPhase;
    PPhase = acqparams.PPhase;
    pPhase = acqparams.pPhase;
    M = acqparams.M; % modulation order
    % search oversampling factor
    searchOS = acqparams.searchOS;
    % search range for Doppler
    f_min = acqparams.f_min;
    f_max = acqparams.f_max;
    % pskacquisition.m precorrects with the midpoint of the Doppler rate range
    frhat = acqparams.fr_min + (acqparams.fr_max-acqparams.fr_min)/2;
    % search range for Doppler rate (CAUTION: 2nd derivative of phase produces a factor 2)
    fr_min = (acqparams.fr_min - frhat)/2;
    fr_max = (acqparams.fr_max - frhat)/2;
    
    if ((fr_min-fr_max) ~= 0)
        mexCoherentMackenthun(phaseEstimatorName, DPhase, PPhase, pPhase, M, f_min, f_max, fr_min, fr_max, T, searchOS);
    else
        if ((f_min-f_max) ~= 0)
            mexCoherentMackenthun(phaseEstimatorName, DPhase, PPhase, pPhase, M, f_min, f_max, T, searchOS);
        else
            mexCoherentMackenthun(phaseEstimatorName, DPhase, PPhase, pPhase, M);
        end
    end
end

end