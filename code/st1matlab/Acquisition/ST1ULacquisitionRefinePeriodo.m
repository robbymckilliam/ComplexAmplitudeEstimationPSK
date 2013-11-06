% Least squares estimator of frequ., frequ. rate and complex amplitude.
%-------------------------------------------------------------
% ST1ULacquisitionRefinePeriodo.m 
%
% [fhat,frhat,amplhat] = ST1ULacquisitionRefinePeriodo(Rx,TxHat,tSamp,frhatcoarse,AcqParams)
%
% This function implements the least squares etimator of the complex
% amplitude, Doppler frequency and rate of the received signal. The
% estimator is related to the periodogram. First, coarse Doppler frequency
% and rate  estimates are obtained by computing the periodogram at the
% discrete frequencies of an FFT, if required for multiple hypothesised
% Doppler rates. The coarse estimates then serve as initialisation for the
% second estimator stage, which performs a joint hill climb on the
% objective function of both variables using Newton's method. Similar to
% [1], Newton's method is applied to the logarithm of the periodogram
% rather than the periodogram itself. This approach eliminates the need for
% zero padding to increase the resolution of the FFT-based initial
% estimate. Finally a complex gain estimate is obtained.
%
%   [1] B. Quinn et al., "Maximizing the Periodogram", in Proc. IEEE Global
%       Global Telecommunications Conference, Nov. 2008 
% 
% Inputs:
%   Rx          - Vector containing samples of received signal
%   TxHat       - Vector containing samples of estimated transmit signal
%                 (typically this signal is rebuild with the soft output of
%                 the decoder)
%   tSamp       - Vector of sample times in seconds
%   AcqParams   - Struct holding aquisition parameters
%       .T          - symbol period
%       .OS         - oversampling factor
%       .fr_min_refine     - lower limit for Doppler rate search range [Hz/s]
%       .fr_max_refine     - upper limit for Doppler rate search range [Hz/s]
%
% Outputs:
%   fhat        - Frequency estimate [Hz]
%   frhat       - Estimate of the Doppler rate [Hz/s]
%   amplhat     - Complex gain estimate
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

function [fhat,frhat,amplhat] = ST1ULacquisitionRefinePeriodo(Rx,TxHat,tSamp,frhatcoarse,AcqParams)

y = Rx .* conj(TxHat);
sumRxHat2 = sum(abs(TxHat).^2);
Ts = AcqParams.T/AcqParams.OS;
N = length(Rx);

fr_min = AcqParams.fr_min_refine;
fr_max = AcqParams.fr_max_refine;

%% coarse frequency and Doppler rate estimation via FFT-based periodogram
pad = 1; % zero-padding factor for FFT
% Nfft = pad*N;
Nfft = AcqParams.Nfft;
if isempty(Nfft)
    Nfft = 2^nextpow2(pad*N);
elseif (Nfft < N)
    y       = y(1:Nfft);
    tSamp   = tSamp(1:Nfft);
end

if ((AcqParams.fr_min-AcqParams.fr_max) ~= 0) && ((AcqParams.f_min-AcqParams.f_max) ~= 0)
    
    % Estimate required Doppler rate resolution based on assumption that 4*N^2
    % samples are required to cover the entire range eta\in[-1/4,1/4) (see
    % Robby's thesis).
    eta_reso = 0.5/(pad*N^2)/Ts^2;
    fr_reso = 2*eta_reso; % factor 2 due to Doppler rate = 2nd derivative of (f*t + eta*t^2)
    
    % Select finer resolution. This has not been optimised in regards to the
    % accuracy-complexity trade-off. It should still work when this line is
    % commented out.
    fr_reso = 0.5*fr_reso;
    
    fr = frhatcoarse + linspace(fr_min,fr_max,ceil(abs(fr_max-fr_min)/fr_reso)+1);
    
    % compute periodogram for each Doppler rate hypothesis
    Y = zeros(length(fr),Nfft);
    for idxfr = 1:length(fr)
        % Doppler rate compensated signal
        rxhypo = y .* exp(-2i*pi*(fr(idxfr)/2)*tSamp.^2);
        % compute periodogram
        F = fft(rxhypo,Nfft);
        
        if isfield(AcqParams, 'FFTMask')
            Y(idxfr,:) = F.*fftshift(AcqParams.FFTMask);
        else
            Y(idxfr,:) = F;
        end
    end
    [Imax,idxmax] = max(abs(Y).^2,[],2);
    [~,idxfrmax] = max(Imax);
    frhatcoarse = fr(idxfrmax);
    fhatcoarse  = (idxmax(idxfrmax)-1)/Nfft/Ts;
    
    % refine estimates via Newton's method
    alpha = 0; % alpha=0: max log(periodogram), alpha=1: max standard periodogram
    gamma = Inf(2,1);
    EPSILON = 1e-10;
    MAX_ITER = 15;
    numIter = 0;
    fhat = fhatcoarse;
    etahat = frhatcoarse/2;
    while (any(abs(gamma) > EPSILON) && (numIter < MAX_ITER))
        % evaluate periodogram I and its first two derivatives Idd and Idd
        X = y .* exp(-2i*pi * (fhat*tSamp + etahat*tSamp.^2));
        Y = 1/sumRxHat2*sum(X);
        %     Ydf     =  -2i*pi/sumRxHat2 * sum(tSamp   .*X);
        %     Ydfr    =  -2i*pi/sumRxHat2 * sum(tSamp.^2.*X);
        %     Ydf2    = -4*pi^2/sumRxHat2 * sum(tSamp.^2.*X);
        %     Ydfrdf  = -4*pi^2/sumRxHat2 * sum(tSamp.^3.*X);
        %     Ydfr2   = -4*pi^2/sumRxHat2 * sum(tSamp.^4.*X);
        
        tmp     = tSamp.*X;
        Ydf     =  -2i*pi/sumRxHat2 * sum(tmp);
        tmp     = tmp.*tSamp;
        Ydfr    =  -2i*pi/sumRxHat2 * sum(tmp);
        Ydf2    = -4*pi^2/sumRxHat2 * sum(tmp);
        tmp     = tmp.*tSamp;
        Ydfrdf  = -4*pi^2/sumRxHat2 * sum(tmp);
        tmp     = tmp.*tSamp;
        Ydfr2   = -4*pi^2/sumRxHat2 * sum(tmp);
        
        I       = abs(Y)^2;
        Idf     = 2*real(   Ydf * conj(Y));
        Idfr    = 2*real(  Ydfr * conj(Y));
        Idf2    = 2*real(  Ydf2 * conj(Y)) + 2*abs(Ydf)^2;
        Idfr2   = 2*real( Ydfr2 * conj(Y)) + 2*abs(Ydfr)^2;
        Idfrdf  = 2*real(Ydfrdf * conj(Y)  + Ydf*conj(Ydfr));
        
        G = [Idf^2, Idf*Idfr; Idf*Idfr, Idfr^2];
        H = [Idf2 , Idfrdf  ; Idfrdf  , Idfr2 ];
        
        % Newton iterate with monotonic function kappa(I) = I^alpha.
        % For alpha=1, this is the standard Newton iterate H^(-1)*[Idf; Idfr].
        gamma = ((alpha-1)/I*G + H)\[Idf; Idfr];
        fhat   = fhat   - gamma(1);
        etahat = etahat - gamma(2);
        
        numIter = numIter + 1;
    end
    frhat = 2*etahat;

else
    % assume zero frequency and frequency rate offsets
    frhat = 0;
    fhat  = 0;
end

%% handle negative frequency offsets
if (fhat > 1/(2*Ts))
    % negative frequency offset
    fhat = fhat - 1/Ts;
end

%% complex amplitude estimate (least squares estimator)
% evaluate periodogram for refined estimates
if ~strcmpi(AcqParams.phaseEstimatorName,'disable')
    X = y .* exp(-2i*pi * (fhat*tSamp + frhat/2*tSamp.^2));
    amplhat = 1/sumRxHat2*sum(X);
else
    % assume zero phase offset and unit channel gain
    amplhat = sqrt(AcqParams.SigPowEst);
end