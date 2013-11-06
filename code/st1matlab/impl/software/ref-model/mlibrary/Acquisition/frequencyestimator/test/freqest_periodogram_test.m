% Test script for freqest_periodogram.m.

%% system parameters

% local variables
T = 1/10e3; %symbol period, 10KHz symbol rate
OS = 4;
Ts = T/OS;
taumin = 0;
taumax = 10e-3; %maximum 10 millisconds time of flight 
M = 4; %QPSK
g = @(t) sinc(t/T)/sqrt(T); %transmission pulse
L = 1000; % number of symbols
P = []; % pilot indices
D = 1:L; % data indices

% complex exponential for frequency estimation
fsin = (1+0.1)/(2*T); % frequency shift relative to band centre
signal_to_sine_ratio_dB = 10;
amplsin = sqrt( 1 / (10^(signal_to_sine_ratio_dB/10)*T) );
fr_max = 0; % maximum Doppler rate [Hz/s]

%% channel parameters
a0 = 1;
tau0 = 5.111e-3;
f0 = -2364.15; % frequency offset [Hz]
fr0 = 0; % Doppler rate [Hz/s]

snr_dbs = 3;
sigma2 = 10.^(-snr_dbs/10)/Ts;
N = ceil((T*L+taumax)/Ts);
ts = (0:N-1)*Ts;

%struct with the modulation parameters and estimator parameters    
acqparams.T = T; %symbol period, 10KHz symbol rate
acqparams.OS =  OS;
acqparams.taumin = taumin;
acqparams.taumax = taumax; %maximum 10 millisconds time of flight 
acqparams.M = M; %QPSK
acqparams.pshape = g; %transmission pulse
acqparams.P = P; %pilot indices
acqparams.D = D; %data indices
acqparams.g = g;

acqparams.fsin = fsin;
acqparams.amplsin = amplsin;
acqparams.fr_max = fr_max;

%generate random psk symbols
s = exp(2i*pi*(randi(M,[1,L])-1)/M);
% set the pilots
acqparams.p = s(P);

%generate recieved signal.
received = @(t) a0*transmitted(t, s, acqparams);
% complex exponential
cexp = @(t) amplsin*exp(2i*pi*fsin*t) .* (t>=T) .* (t<L*T);
% apply frequency offset and Doppler rate
rx = (received(ts - tau0) + cexp(ts - tau0)) .* exp(2i*pi * (f0*ts + (fr0/2)*ts.^2));
% add noise
rx = rx + sqrt(sigma2)/sqrt(2)*(randn(1,N) + 1i*randn(1,N));

%run the estimator
fhat = freqest_periodogram(rx,ts,acqparams);
% correct frequency estimate by offset of complex exponential
fhat = fhat - fsin;

% The frequency RMSEE at 3dB for a compl. exp. 10dB below the signal is
% approx. 0.6Hz for a symbol rate of 16.384kHz. Assuming a normal
% distribution for the residual error, 99.73% of the errors are within the
% 3*sigma confidence interval.
tol = 3*0.6/16384;
if(abs(fhat - f0)*T < tol)
    disp('Testing frequency offset estimator. PASS');
%     pass(end+1) = true;
else
    disp('Testing frequency offset estimator. FAIL');
%     pass(end+1) = false;
end