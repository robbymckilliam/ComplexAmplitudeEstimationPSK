% Test script for mexFrequencyEstimator.
%--------------------------------------------------------------------------
% Author:  Andre Pollok, 20/02/2013
% Project: ASRP
%--------------------------------------------------------------------------
% Copyright 2013 : 
% Institute for Telecommunications Research
% University of South Australia
%--------------------------------------------------------------------------

debug_plot = false;

numSamp = 100;
Ts = 1; % sample period

% tolerance for pass tests
tol = 1e-5;

% channel parameter offsets
ftrue = -0.831/(2*Ts);
frtrue = +0.0009;
cgaintrue = 0.75*exp(1i*pi/4);

% estimator settings
Nfft = 128; % FFT size
frmin = -0.001; % lower limit for frequency rate search range in Hz/s
frmax = +0.001; % upper limit for frequency rate search range in Hz/s
searchOS = 1;

% sample times
tauref = 10*Ts;
tSamp = [0:numSamp-1]*Ts - tauref;
% transmitted signal (unit-amplitude tone at DC)
txSamp = ones(1,numSamp);

%% received signal with complex gain and frequency offset
rxSamp = txSamp * cgaintrue .* exp(2i*pi*ftrue*tSamp);

if debug_plot
    figure
    plot(fftshift(abs(fft(rxSamp,Nfft)).^2))
    hold on, grid on
end

fprintf('Testing estimator with frequency offset... ');
% construct and run estimator
mexFrequencyEstimator('freqEst',Nfft,Ts,tauref);
[fhat,cgainhat] = mexFrequencyEstimator('freqEst',rxSamp,txSamp);
% [fhat,cgainhat] = mexFrequencyEstimator('freqEst',rxSamp); % implicitly assumes unit-amplitude tone at DC
pass = (abs(ftrue - fhat) < tol) && (abs(cgaintrue - cgainhat) < tol);
if pass
    pass = 'PASS';
else
    pass = 'FAIL';
end
fprintf([pass,'\n'])

%% received signal with complex gain, frequency and frequency rate offset
rxSamp = txSamp * cgaintrue .* exp(2i*pi * (ftrue*tSamp + (frtrue/2)*tSamp.^2));

if debug_plot
    plot(fftshift(abs(fft(rxSamp,Nfft)).^2), 'r-')
end

fprintf('Testing estimator with frequency and frequency rate offsets... ');
% construct and run estimator
mexFrequencyEstimator('freqEst',Nfft,Ts,tauref,frmin,frmax,searchOS);
[fhat,cgainhat,frhat] = mexFrequencyEstimator('freqEst',rxSamp,txSamp);
% [fhat,cgainhat,frhat] = mexFrequencyEstimator('freqEst',rxSamp); % implicitly assumes unit-amplitude tone at DC
pass = (abs(ftrue - fhat) < tol) && (abs(frtrue - frhat) < tol) && (abs(cgaintrue - cgainhat) < tol);
if pass
    pass = 'PASS';
else
    pass = 'FAIL';
end
fprintf([pass,'\n'])