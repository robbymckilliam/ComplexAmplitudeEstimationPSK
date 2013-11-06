clear all;

addpath ../bin
addpath ../..
    
%%%%%%% Phase only test 

%local variables
T = 1/10e3; %symbol period, 10KHz symbol rate
M = 4; %QPSK
L = 200; %number of symbols
absP = 50; %number of pilots
P = 1:absP; %pilot indices
D = absP+1:L; %data indices

s = exp(2*pi*1i*(randi(M,[1,L])-1)/M); %qpsk symbols
p = s(P); %get the pilot symbols

a0 = 1.0 + 0.0i; %true complex gain

%create our estimator
mexCoherentMackenthun('CoherenMackenthun', D, P, p, M);

%create params struct for the matlab version
params.M =  M;
params.D = D;
params.P = P;
params.p = p;
params.L =  L;

snrsdb = 0:10:80;
snrs = 10.^(snrsdb/10);
vars = abs(a0)^2/2./snrs; %variance of real and imaginary parts.

iters = 10;
pass = true;
TOL = 1e-6;

for snritr = 1:length(snrs)
    
    errsa = zeros(1,iters);
    errsf = zeros(1,iters);
    for itr = 1:iters
        w = sqrt(vars(snritr))*(randn(1,L) + 1i*randn(1,L)); %noise
        y = a0*s + w; %received signal
        
        ahat = mexCoherentMackenthun('CoherenMackenthun', y);
        ahattest =  coherentmackenthun(y, params);
        %abs(ahattest - ahat)
        
        pass = pass && abs(ahattest - ahat) < TOL;
    end

end


if(pass) disp('Testing Coherent Mackenthun ... PASS');
else disp('Testing Coherent Mackenthun ... FAIL');
end
    

    
    
%%%%%%% With Doppler only test    
    
%local variables
T = 1/10e3; %symbol period, 10KHz symbol rate
M = 4; %QPSK
L = 100; %number of symbols
absP = 20; %number of pilots
P = 1:absP; %pilot indices
D = absP+1:L; %data indices

s = exp(2*pi*1i*(randi(M,[1,L])-1)/M); %qpsk symbols
p = s(P); %get the pilot symbols

a0 = 1.0 + 0.0i; %true complex gain
f0 = 4.7;

fmin = -5;
fmax = 5;

%create our estimator
mexCoherentMackenthun('CoherenMackenthunWithDoppler', D, P, p, M, ...
                      fmin, fmax, T);

%create params struct for the matlab version
params.M =  M;
params.D = D;
params.P = P;
params.p = p;
params.L =  L;
params.T = T;
params.f_min = fmin;
params.f_max = fmax;

snrsdb = 0:20:80;
snrs = 10.^(snrsdb/10);
vars = abs(a0)^2/2./snrs; %variance of real and imaginary parts.

iters = 5;
ts = T*(1:L);

pass = true;
TOL = 1e-6;

for snritr = 1:length(snrs)
    
    errsa = zeros(1,iters);
    errsf = zeros(1,iters);
    for itr = 1:iters
        w = sqrt(vars(snritr))*(randn(1,L) + 1i*randn(1,L)); %noise
        y = a0*s.*exp(2i*pi*ts*f0) + w; %received signal
        
        [ahat, fhat] = ...
            mexCoherentMackenthun('CoherenMackenthunWithDoppler', y);
        [ahattest, fhattest] =  coherentmackenthunwithonlydoppler(y, ...
                                                          params);
        
        
        pass = pass && abs(ahattest - ahat) < TOL;
        pass = pass && abs(fhattest - fhat) < TOL;
        %if(abs(fhattest - fhat) > TOL)
        %    disp([ snrsdb(snritr), abs(fhattest - f0), abs(fhat - ...
        %                                                   f0) ]);
        %end
    end

end


if(pass) disp('Testing Coherent Mackenthun with only Doppler ... PASS');
else disp('Testing Coherent Mackenthun with only Doppler ... FAIL');
end


 %%%%%%% With Doppler only and noncontiguous symbols test    
     
 %local variables
 T = 1/10e3; %symbol period, 10KHz symbol rate
 M = 4; %QPSK
 L = 100; %number of symbols
 P = 1:20; %pilot indices
 D = 50:80; %data indices
 
 s = exp(2*pi*1i*(randi(M,[1,L])-1)/M); %qpsk symbols
 p = s(P); %get the pilot symbols
 
 a0 = 1.0 + 0.0i; %true complex gain
 f0 = 4.7;
 
 fmin = -5;
 fmax = 5;
 
 %create our estimator
 mexCoherentMackenthun('CoherenMackenthunWithDoppler', D, P, p, M, ...
                       fmin, fmax, T);
 
 %create params struct for the matlab version
 params.M =  M;
 params.D = D;
 params.P = P;
 params.p = p;
 params.L =  L;
 params.T = T;
 params.f_min = fmin;
 params.f_max = fmax;
 
 snrsdb = 0:20:80;
 snrs = 10.^(snrsdb/10);
 vars = abs(a0)^2/2./snrs; %variance of real and imaginary parts.
 
 iters = 5;
 ts = T*(1:L);
 
 pass = true;
 TOL = 1e-6;
 
 for snritr = 1:length(snrs)
     
     errsa = zeros(1,iters);
     errsf = zeros(1,iters);
     for itr = 1:iters
         w = sqrt(vars(snritr))*(randn(1,L) + 1i*randn(1,L)); %noise
         y = a0*s.*exp(2i*pi*ts*f0) + w; %received signal
         
         [ahat, fhat] = ...
             mexCoherentMackenthun('CoherenMackenthunWithDoppler', y);
         [ahattest, fhattest] =  coherentmackenthunwithonlydoppler(y, ...
                                                           params);
         
         
         pass = pass && abs(ahattest - ahat) < TOL;
         pass = pass && abs(fhattest - fhat) < TOL;
         if(abs(fhattest - fhat) > TOL)
             disp([ snrsdb(snritr), abs(fhattest - f0), abs(fhat - ...
                                                           f0) ]);
         end
     end
 
 end
 
 
 if(pass) disp(['Testing Coherent Mackenthun with Doppler and ' ...
                'noncontiguous symbols... PASS']);
 else disp(['Testing Coherent Mackenthun with Doppler and noncontiguous ' ...
            'symbols ... FAIL']);
end


%%%%%% Doppler and Doppler rate test
    
%local variables
T = 1/10e3; %symbol period, 10KHz symbol rate
M = 4; %QPSK
L = 200; %number of symbols
absP = 100; %number of pilots
P = 1:absP; %pilot indices
D = absP+1:L; %data indices

s = exp(2*pi*1i*(randi(M,[1,L])-1)/M); %qpsk symbols
p = s(P); %get the pilot symbols

a0 = cos(pi/4) + i*sin(pi/4); %true complex gain
f0 = 4.7;
fr0 = -11.1;

fmin = -5;
fmax = 5;
frmin = -30;
frmax = 30;

%create our estimator
mexCoherentMackenthun('CoherenMackenthunWithRate', D, P, p, M, ...
                      fmin, fmax, frmin, frmax, T);

%create params struct for the matlab version
params.M =  M;
params.D = D;
params.P = P;
params.p = p;
params.L =  L;
params.T = T;
params.f_min = fmin;
params.f_max = fmax;
params.fr_min = frmin;
params.fr_max = frmax;

snrsdb = 0:20:80;
snrs = 10.^(snrsdb/10);
vars = abs(a0)^2/2./snrs; %variance of real and imaginary parts.

iters = 5;
ts = T*(1:L);

pass = true;
TOL = 1e-6;

for snritr = 1:length(snrs)
    
    errsa = zeros(1,iters);
    errsf = zeros(1,iters);
    for itr = 1:iters
        w = sqrt(vars(snritr))*(randn(1,L) + 1i*randn(1,L)); %noise
        y = a0*s.*exp(2i*pi*ts*f0 + 2i*pi*ts.^2*fr0) + w; %received signal
        
        [ahat, fhat,frhat] = ...
            mexCoherentMackenthun('CoherenMackenthunWithRate', y);
        [ahattest, fhattest, frhattest] = ...
            coherentmackenthunwithdoppler(y, params);
        
        %if(abs(ahattest - ahat) > TOL) 
        %    disp([abs(ahat - a0), abs(ahattest - a0)]);
        %end
        %disp(angle(ahat*conj(ahattest)));
        
        pass = pass && abs(ahattest - ahat) < TOL;
        pass = pass && abs(fhattest - fhat) < TOL;
        pass = pass && abs(frhattest - frhat) < TOL;
        if(~pass)
            disp([ snrsdb(snritr), abs(fhattest - f0), abs(fhat - ...
                                                           f0) ]);
            disp([ snrsdb(snritr), abs(frhattest - fr0), abs(frhat - ...
                                                           fr0) ]);
        end
    end

end


if(pass) disp('Testing Coherent Mackenthun with Doppler Rate ... PASS');
else disp('Testing Coherent Mackenthun with Doppler Rate ... FAIL');
end


%%%%%% Doppler and Doppler rate test with noncontiguous symbols
    
 %local variables
 T = 1/10e3; %symbol period, 10KHz symbol rate
 M = 4; %QPSK
 L = 100; %number of symbols
 P = 1:20; %pilot indices
 D = 50:80; %data indices

s = exp(2*pi*1i*(randi(M,[1,L])-1)/M); %qpsk symbols
p = s(P); %get the pilot symbols

a0 = cos(pi/4) + i*sin(pi/4); %true complex gain
f0 = 4.7;
fr0 = -11.1;

fmin = -5;
fmax = 5;
frmin = -30;
frmax = 30;

%create our estimator
mexCoherentMackenthun('CoherenMackenthunWithRate', D, P, p, M, ...
                      fmin, fmax, frmin, frmax, T);

%create params struct for the matlab version
params.M =  M;
params.D = D;
params.P = P;
params.p = p;
params.L =  L;
params.T = T;
params.f_min = fmin;
params.f_max = fmax;
params.fr_min = frmin;
params.fr_max = frmax;

snrsdb = 0:20:80;
snrs = 10.^(snrsdb/10);
vars = abs(a0)^2/2./snrs; %variance of real and imaginary parts.

iters = 5;
ts = T*(1:L);

pass = true;
TOL = 1e-6;

for snritr = 1:length(snrs)
    
    errsa = zeros(1,iters);
    errsf = zeros(1,iters);
    for itr = 1:iters
        w = sqrt(vars(snritr))*(randn(1,L) + 1i*randn(1,L)); %noise
        y = a0*s.*exp(2i*pi*ts*f0 + 2i*pi*ts.^2*fr0) + w; %received signal
        
        [ahat, fhat,frhat] = ...
            mexCoherentMackenthun('CoherenMackenthunWithRate', y);
        [ahattest, fhattest, frhattest] = ...
            coherentmackenthunwithdoppler(y, params);
        
        %if(abs(ahattest - ahat) > TOL) 
        %    disp([abs(ahat - a0), abs(ahattest - a0)]);
        %end
        %disp(angle(ahat*conj(ahattest)));
        
        pass = pass && abs(ahattest - ahat) < TOL;
        pass = pass && abs(fhattest - fhat) < TOL;
        pass = pass && abs(frhattest - frhat) < TOL;
        if(~pass)
            disp([ snrsdb(snritr), abs(fhattest - f0), abs(fhat - ...
                                                           f0) ]);
            disp([ snrsdb(snritr), abs(frhattest - f0), abs(frhat - ...
                                                           f0) ]);
         end
    end

end


if(pass) disp(['Testing Coherent Mackenthun with Doppler Rate with ' ...
               'noncontiguous symbols ... PASS']);
else disp(['Testing Coherent Mackenthun with Doppler Rate with ' ...
           'noncontiguous symbols ... FAIL']);
end


%%%%%% Test correctly deleted

%local variables
T = 1/10e3; %symbol period, 10KHz symbol rate
M = 4; %QPSK
L = 200; %number of symbols
absP = 100; %number of pilots
P = 1:absP; %pilot indices
D = absP+1:L; %data indices

s = exp(2*pi*1i*(randi(M,[1,L])-1)/M); %qpsk symbols
p = s(P); %get the pilot symbols

a0 = cos(pi/4) + i*sin(pi/4); %true complex gain
f0 = 4.7;
fr0 = -11.1;

fmin = -5;
fmax = 5;
frmin = -30;
frmax = 30;

itrs = 20000; %reconstruct 20000 times

for itr = 1:itrs
    
    mexCoherentMackenthun('CoherenMackenthun', D, P, p, M);
    
    mexCoherentMackenthun('CoherenMackenthunWithDoppler', D, P, p, M, ...
                    fmin, fmax, T);

    mexCoherentMackenthun('CoherenMackenthunWithRate', D, P, p, M, ...
                     fmin, fmax, frmin, frmax, T);
end

%if we get here, we have passed
disp('Testing construction and deletion ... PASS');
