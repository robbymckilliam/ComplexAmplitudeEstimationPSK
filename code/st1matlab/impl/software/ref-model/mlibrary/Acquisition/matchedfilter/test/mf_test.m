addpath ..
%addpath ../../timeoffsetestimator/test
%addpath ../../timeoffsetestimator

%local variables
T = 1/10e3; %symbol period, 10KHz symbol rate
OS = 4;
Ts = T/OS;
taumin = 0;
taumax = 15e-3; %maximum 10 millisconds time of flight 
p = 1;
q = OS;
c = 4; 
a = 4;
Delta = Ts/c;
taugrid = taumin:Delta:taumax;
K = length(taugrid);
L = 1000;
N = ceil((T*L+taumax)/Ts);
g = @(t) sinc(t/T)/sqrt(T); %transmission pulse
%figure(1); ts = -5:0.1:5; plot(ts, g(ts));
s = exp(2*pi*1i*(randi(4,[1,L])-1)/4); %qpsk symbols
absP = 10; %number of pilots
P = 1:absP; %pilot indices
D = absP+1:L; %data indices
maxD = max(D); %maximum data index
minD = min(D); %maximum data index

%struct with the modulation parameters and esitmator parameters    
params.T = T ; %symbol period
params.p = p ; 
params.q = q ;
params.OS = q;
params.Ts = Ts; %sample period
params.taumin = taumin;
params.taumax = taumax ;
params.c = c ; 
params.a = a ;
params.Delta = Delta;
params.taugrid = taugrid;
params.K = K;
params.L = L;
params.N = N;
params.g = g; %transmission pulse
%figure(1); ts = -5:0.1:5; plot(ts, g(ts));
params.s = s; %mpsk symbols
params.pilots = s(P); %mpsk pilots
params.absP = absP; %number of pilots
params.P = P; %pilot indices
params.D = D; %data indices
params.maxD = maxD; %maximum data index
params.minD = minD; %minimum data index
params.maxP = max(P); %maximum data index
params.minP = min(P); %minimum data index
params.pshape = g;

%true timeoffset and phase
tau0 = 5.111e-3; 
%a0 = exp(2*pi*1i*rand);
a0 = 1;

%signal to noise ration and variance
SNRdB = 20;
SNR = 10^(SNRdB/10);
sigma2 = Ts/T/SNR;

%received signal
recieved = @(t) a0*transmitted(t - tau0, s, params);
r = recieved(Ts*(1:N)) + sqrt(sigma2)/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
%figure(2); ts = -30*T:Delta:Ts*N+30*T; plot(ts, abs(transmitted(ts, s, params))); hold on;
%plot(ts, abs(recieved(ts)), 'r');
%plot(Ts*(1:N), abs(r), 'g.');
%title('transmitted and recieved signals'); xlabel('time t'); ...
%    ylabel('signal amplitude');
%csvwrite('xreal',real(x(ts)))
%csvwrite('ximag',imag(x(ts)))
%csvwrite('rreal',real(recieved(ts)))
%csvwrite('rimag',imag(recieved(ts)))

%add the recieved signal to the params struct
params.r = r;

mfout = mf(r, tau0, params);

%figure(6); plot(mfout - s, '.'); hold on;
%plot(s,'.r');
%plot(mfout, '.g');

 if( sum((mfout - s).^2)/length(s) < 0.05 )
     disp('Testing matched filter.  PASS')
 else
     disp('Testing matched filter.  FAIL')
 end
