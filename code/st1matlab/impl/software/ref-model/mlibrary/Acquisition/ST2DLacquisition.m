function [mfoutput, fhat, frhat, cgainhat, tauhat, nvarhat] = ...
    ST2DLacquisition(rx, acqparams)
% This function provides channel parameter offset estimates for ST2 DL PSK 
% signals. The function operates on an oversampled received signal and
% produces estimates for carrier frequency offset, frequecy rate, timing
% offset, complex signal amplitude and noise variance. Also returns a
% matched filtered output at symbol rate.
%
% COMMENT: the time reference point (t=0) for the returned estimates of the
%          polynomial phase coefficients is the midpoint of the transmitted
%          burst, i.e. (tauhat + L/2*T), where L is the number of symbols
%          in the burst.
%
%   [mfoutput, fhat, frhat, cgainhat, tauhat, nvarhat] = ST1ULacquisition(rx,acqparams)
%
% INPUTS:
%   rx              - sampled received signal
%   acqparam        - struct holding aquisition parameters
%       .T              - symbol period
%       .OS             - oversampling factor
%       .p              - pilot symbols
%       .M              - M for M-PSK, i.e. size of constellation
%       .fsin           - frequency offset of the sinusoid
%       .pshape         - function representing transmission pulse shape
%       .pshapeName     - string identifying the pulse shape
%       .rolloff        - roll-off factor of root raised cosine pulse
%       .Tpshape        - duration of the transmit pulse shape in seconds
%       .timingEstimatorName
%                       - string identifying the timing estimator
%       .phaseEstimatorName
%                       - string identifying the phase estimator
%       .fr_max         - maximum Doppler rate
%       .fr_min         - minimum Doppler rate
%       .freqEstimatorName
%                       - string identifying the initial frequency
%                         estimator (periodogram)
%       .Nfft           - FFT size for initial frequency estimator
%
% OUTPUTS:
%   mfoutput        - matched filtered signal at symbol rate
%   fhat            - frequency estimate
%   frhat           - frequency rate estimate
%   tauhat          - timing estimate
%   cgainhat        - complex signal amplitude estimate (gain and phase)
%   nvarhat         - Estimate of the noise variance per symbol and signal
%                     dimension. Assumes that the real and imaginary parts
%                     of the noise have the same variance. The Es/No can be
%                     computed as abs(cgainhat)^2/(2*nvarhat).
%
% (c) 2012 ITR-UniSA
% Authors: Andre Pollok
% Created: 2 August 2012 (based on pskacquisition.m)

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
    acqparams.p = acqparams.p * exp(1i*pi/4);
end

%% Estimate frequency offset from periodogram
T = acqparams.T; % symbol period
OS = acqparams.OS; % oversampling factor
Ts = T/OS; % sample period
N = length(rx);
ts = (0:N-1)*Ts; % sample times

% precorrection with the midpoint of the Doppler rate range
frhat = acqparams.fr_min + (acqparams.fr_max-acqparams.fr_min)/2;
rxcorr = rx .* exp(-2i*pi*(frhat/2)*ts.^2);

% estimate frequency offset (discard samples before minimum time delay)
freqEstimatorName = acqparams.freqEstimatorName;
Nfft = acqparams.Nfft;  % FFT size for periodogram
idxStart = floor(acqparams.taumin/Ts)+1;
if (Nfft >= N-idxStart+1)
    [fhat] = mexFrequencyEstimator(freqEstimatorName,rxcorr(idxStart:end));
else
    [fhat] = mexFrequencyEstimator(freqEstimatorName,rxcorr(idxStart+[0:Nfft-1]));
end

% correct frequency estimate by offset of complex exponential
fsin = acqparams.fsin;
fhat = fhat - fsin;

% frequency correction
rxcorr = rxcorr .* exp(-2i*pi*fhat*ts);

%% Estimate time offset from frequency corrected samples

timingEstimatorName = acqparams.timingEstimatorName;
tauhatMF = mexTimingEstimator(timingEstimatorName,rxcorr);
% adjust time offset estimate by (OS-1) samples (timing estimator and
% matched filter assume the first sample at t=Ts and the first symbol at 
% t=T=OS*Ts rather than at t=0)
tauhat = tauhatMF + (OS-1)*Ts;

%% Filter out tone and run the matched filter (samples to symbols)

if (fsin == 0)
    % remove DC component from frequency-corrected signal
    [rxcorr] = dcblockST2DL(rxcorr);
end

%mfoutput = mf(rxcorr, tauhatMF, acqparams);
mfoutput = mf_coder(rxcorr, tauhatMF, acqparams.rolloff, acqparams.OS, acqparams.T, acqparams.P, acqparams.D, acqparams.Tpshape/2/T);

%% Estimate phase and gain from the matched filtered signal

phaseEstimatorName = acqparams.phaseEstimatorName;
[cgainhatfine, fhatfine, frhatfine] = mexCoherentMackenthun(phaseEstimatorName, mfoutput);
frhatfine = 2*frhatfine; % CAUTION: 2nd derivative of phase produces a factor 2

%% change time point of reference to midpoint of burst

% change time reference point of coarse estimates by tauhat (signal has
% been time-corrected) plus L/2*T
L = length(mfoutput);
[phihat,fhat,frhat] = polynomialPhaseTransform(0,fhat,frhat,tauhat+L/2*T);

% change time reference point of fine estimates by T (CoherenMackenthunWithDoppler
% assumes the first symbol at t=T rather than at t=0) plus L/2*T
[phihatfine,fhatfine,frhatfine] = polynomialPhaseTransform(angle(cgainhatfine),fhatfine,frhatfine,T+L/2*T);
cgainhatfine = abs(cgainhatfine) * exp(1i*phihatfine);

% combine fine and coarse estimates
frhat = frhat + frhatfine;
fhat  = fhat  + fhatfine;
cgainhat = cgainhatfine * exp(1i*phihat);

%% correct with refined frequency estimate

if (fsin ~= 0)
    % aplly fine estimates (coarse estimates have already been corrected for)
    tsym = [0:L-1]*T - L/2*T;
    mfoutput = mfoutput/cgainhatfine .* exp(-2i*pi * (fhatfine*tsym + (frhatfine/2)*tsym.^2));
else
    % DC tone has to be filtered out
    % -> operate at sample level and re-run matched filter
    
    % correct polynomial phase signal (time reference point at tauhat+L/2*T)
    tauref = tauhat + L/2*T;
    rxcorr = rx .* exp(-1i*angle(cgainhat) - 2i*pi * (fhat*(ts-tauref) + (frhat/2)*(ts-tauref).^2));
    
    % remove DC component from frequency-corrected signal
    [rxcorr] = dcblockST2DL(rxcorr);
    
    % re-run matched filter
    mfoutput = mf(rxcorr, tauhatMF, acqparams);
    % normalize MF output
    mfoutput = mfoutput / abs(cgainhat);
end

%% Estimate noise variance

% estimate Es/No
EsNohat = EsNoest(mfoutput,acqparams);
%EsNohat = EsNoest_coder(mfoutput,acqparams.M,acqparams.D,acqparams.P,acqparams.p);
% noise variance per symbol and signal dimension
nvarhat = abs(cgainhat)^2/(2*EsNohat);

%% correct MF output and complex gain estimate to match input constellation
if (acqparams.M == 4)
    % compensate 45 degree mismatch
    cgainhat = cgainhat * exp(+1i*pi/4);
    mfoutput = mfoutput * exp(-1i*pi/4);
end

%% adjust time offset estimate by leading matched filter tail
tauhat = tauhat - acqparams.Tpshape/2;

end