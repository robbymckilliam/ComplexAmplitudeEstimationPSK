addpath ../bin
addpath ..
addpath ./timeoffsettests/oldtimingestimator


%local variables
T = 1/10e3; %symbol period, 10KHz symbol rate
p = 1;
q = 4;
OS = q/p;
Ts = T/OS;
taumin = 0;
taumax = 30e-3; %maximum 10 millisconds time of flight 
M = 4; %QPSK

% root raised cosine
rolloff = 0.5;
Nsym = 10; % filter order in symbols
duration = Nsym*T; % duration of transmit pulse shape
g = @(t) rrcpulse(t,rolloff,T,OS) / sqrt(Ts/OS) .* (abs(t) <= duration/2);

L = 100; %number of symbols
absP = 50; %number of pilots
P = 1:absP; %pilot indices
D = absP+1:L; %data indices

a0 = 1;

%symbols
tau0 = 5.111e-3;
sigma2 = 0.00000001/Ts;
N = ceil((T*L+taumax)/Ts);

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
acqparams.Tpshape = duration;

%generate random psk symbols
s = exp(2*pi*1i*(randi(M,[1,L])-1)/M);
%set the pilots
pilots = s(P);
acqparams.p = pilots;
acqparams.pilots = pilots;

%generate recieved signal.
recieved = @(t) a0*transmitted(t - tau0, s, acqparams);
rx = recieved(Ts*(1:N)) + sqrt(sigma2)/sqrt(2)*(randn(1,N) + ...
                                                  1i*randn(1,N));

%figure(3); plot(Ts*(1:N), rx); hold on;
%plot(Ts*(1:N), recieved(Ts*(1:N)), 'r');


%%%%%%%%%%%%% Construction and deletion tests %%%%%%%%%%%%%%%%
iters = 50;

for n = 1:iters

    mexTimingEstimator('ST1TimingEstimator',P,D,pilots,T,Ts,p,q,taumin,taumax,rolloff,duration);
    mexTimingEstimator('ST2TimingEstimator',P,D,pilots,T,Ts,taumin,taumax,rolloff,duration);

end

disp('Construction and deletion test. ... PASS');
    

%%%%%%%%%%%%%% Time offset estimator tests %%%%%%%%%%%%%%%%%%


    % construct pulse shape and estimator
    mexTimingEstimator('TimingEstimator',P,D,pilots,T,Ts,p,q,taumin,taumax,rolloff,duration);
    % construct pulse shape and estimator
    mexTimingEstimator('TimingEstimatortol',P,D,pilots,T,Ts,p,q,taumin,taumax,rolloff,duration,1e-8);
    % run the estimator
    %tau0
    tauhat = mexTimingEstimator('ST1TimingEstimator',rx) ;
    tauhattol = mexTimingEstimator('TimingEstimatortol',rx);
    tauhatst2 = mexTimingEstimator('ST2TimingEstimator',rx) ;
    
    if(abs(tauhat - tau0) < 1e-6 && abs(tauhattol - tau0) && ...
       abs(tauhatst2 - tau0) < 1e-6)
        disp('Estimate test ... PASS');
    else
        disp('Estimate test ... FAIL');
    end
    
    %%%%%%%%%%%%%%%%% Same as old version %%%%%%%%%%%%%%%%%%%
    
    % construct old estimator
    mexOldTimeEstimator('RootRaisedCosine',rolloff,T,OS,duration);
    mexOldTimeEstimator('TimingEstimator','RootRaisedCosine',D,P,pilots,taumin,taumax);
    % run old estimator
    [tauhatold,fval] = mexOldTimeEstimator('TimingEstimator',rx);
    
    %tauhatold
    %tauhat
    %tau0
    
    if( abs(tauhat - tauhatold) < 1e-6 )
        disp('Same as old estimator ... PASS');
    else
        disp('Same as old estimator ... FAIL');
    end
    
    
    %%%%%%%%%%%%%%%%% ST1 Doesn't crash with short recieve signal %%%%%%%%%%%%%%%%%%%
    
    shortrx = [1,2,3,4];
    
     tauhat = mexTimingEstimator('ST1TimingEstimator',shortrx) ;
     
    disp('ST1 does not crash with short recieve signal ... PASS');
    %disp('ST1 does not crash with short recieve signal ... FAIL');
    
    %%%%%%%%%%%%%%%%% ST2 Doesn't crash with short recieve signal %%%%%%%%%%%%%%%%%%%
    
    shortrx = [1,2,3,4];
    
     tauhat = mexTimingEstimator('ST2TimingEstimator',shortrx) ;

     disp('ST2 does not crash with short recieve signal ... PASS');
     %disp('ST2 does not crash with short recieve signal ... FAIL');