function mfout = mf(rx, tauhat, acqparam)
%Run a matched filter on the received samples, returns recieved
%symbols.
%    
%    mfout = mf(rx, tauhat, acqparams)
%
% INPUTS:
%   rx              - sampled received signal
%   tauhat        - time offset estimate
%   acqparam        - struct holding aquisition parameters
%       .T              - symbol period
%       .L              - number of symbols (assumed consequtive)
%       .OS             - oversampling factor
%       .pshape         - function representing transmission pulse shape
%
% OUTPUTS:
%   mfout        - matched filtered signal at symbol rate
    
% (c) 2012 ITR-UniSA
% Authors: Robby McKilliam
% Created: 4 April 2012
    
    g = acqparam.pshape;
    p = 1;
    q = acqparam.OS;
    T = acqparam.T;
    Ts = p*T/q;
    mini = min(min(acqparam.P),min(acqparam.D)); 
    maxi = max(max(acqparam.P),max(acqparam.D));
    N = length(rx);
    L = maxi;
    
    ns = p*q*mini - p*N : p*q*maxi - p;
    lvec = conj( g(-ns*Ts/p - tauhat) );
    zvec = zeros(1,p*N - p + 1); zvec(1:p:p*N-p+1) = rx; %zerofilled r
    h = conv_fft(lvec, zvec);

    %plot(real(h));
    %length(h)
    %length(q*(mini:maxi))
    
    mfout = h(q*(mini:maxi) - q*mini + 1)*Ts;
    
end