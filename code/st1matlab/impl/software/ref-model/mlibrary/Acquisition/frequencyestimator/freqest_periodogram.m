% Frequency estimation from periodogram.
%-------------------------------------------------------------
% freqest_periodogram.m 
%
% [fhat,cgainhat] = freqest_periodogram(RxSamp,tSamp,AcqParams)
%
% This function etimates the frequency of a noisy complex exponential by
% maximising the periodogram. First, a coarse estimate is obtained by
% computing the periodogram at the discrete frequencies of an FFT. This
% coarse estimate then serves as the initialisation for the second
% estimator stage, which climbs the objective function using Newton's
% method. As proposed by Quinn et al. [1], Newton's method is applied to
% the logarithm of the periodogram rather than the periodogram itself. This
% approach eliminates the need for zero padding to increase the resolution
% of the FFT-based initial estimate.
%
%   [1] B. Quinn et al., "Maximizing the Periodogram", in Proc. IEEE Global
%       Global Telecommunications Conference, Nov. 2008 
% 
% Inputs:
%   RxSamp      - Vector containing samples of received signal
%   tSamp       - Vector of sample times in seconds
%   AcqParams   - Struct holding aquisition parameters
%       .T          - symbol period
%       .OS         - oversampling factor
%       .taumin     - minimum time offset
%       .Nfft       - FFT length. If set to [] the length of RxSamp rounded
%                     up to the next higher power of 2 will be used.
%
% Outputs:
%   fhat        - Frequency estimate [Hz]
%   cgainhat    - complex gain estimate
%
% Comments:
% 
%-------------------------------------------------------------
% Author: A. Pollok, Institute for Telecommunications Research
% Project: ASRP
%-------------------------------------------------------------
% Copyright 2012 
% Institute for Telecommunications Research
% University of South Australia
%-------------------------------------------------------------

function [fhat,cgainhat] = freqest_periodogram(RxSamp,tSamp,AcqParams)

%% system parameters
T = AcqParams.T; % symbol period
OS = AcqParams.OS; % oversampling factor
Ts = T/OS; % sample period

% channel parameters
taumin = AcqParams.taumin;

% discard samples before minimum time delay
idxStart = floor(taumin/Ts)+1;
RxSamp = RxSamp(idxStart:end);
tSamp  =  tSamp(idxStart:end);

% set FFT length for periodogram
N = length(RxSamp); % number of samples
Nfft = AcqParams.Nfft;
if isempty(Nfft)
    pad = 1; % zero-padding factor for FFT
    Nfft = 2^nextpow2(pad*N);
elseif (Nfft < N)
    RxSamp = RxSamp(1:Nfft);
    tSamp  = tSamp(1:Nfft);
end

%% coarse frequency estimation via periodogram
I = fft(RxSamp,Nfft);

if isfield(AcqParams, 'FFTMask')
  I = I.*fftshift(AcqParams.FFTMask);
end

[~,idxmax] = max(abs(I).^2);
fhatcoarse = (idxmax-1)/Nfft/Ts;

%% refine estimate via Newton's method
alpha = 0; % alpha=0: max log(periodogram), alpha=1: max standard periodogram
gamma = Inf;
EPSILON = 1e-10;
MAX_ITER = 15;
numIter = 0;
fhat = fhatcoarse;
while ((abs(gamma) > EPSILON) && (numIter < MAX_ITER))
    % evaluate periodogram I and its first two derivatives Idd and Idd
    X = RxSamp.*exp(-2i*pi*fhat*tSamp);
    Y = 1/Nfft*sum(X);
    Yd =  -2i*pi/Nfft  * sum(tSamp   .*X);
    Ydd = -4*pi^2/Nfft * sum(tSamp.^2.*X);
    I = abs(Y)^2;
    Id =  2*real(Yd*conj(Y));
    Idd = 2*real(Ydd*conj(Y)) + 2*abs(Yd)^2;
    % Newton iterate
    gamma = Id/(Idd + (alpha-1)*Id^2/I);
    fhat = fhat - gamma;
    numIter = numIter + 1;
end

%% handle negative frequency offsets
if (fhat > 1/(2*Ts))
    % negative frequency offset
    fhat = fhat - 1/Ts;
end

%% estimate complex amplitude (least squares estimator)
if (nargout > 1)
    cgainhat = Y;
end

end