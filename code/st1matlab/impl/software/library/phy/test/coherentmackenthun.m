function [ahat, Qhat] = coherentmackenthun(rx, acqparams)
% This function provides channel parameter offset estimates for PSK 
% signals. The function operates on an oversampled received
% signal and produces estimates for carrier frequency offset,
% frequecy rate timing offset, complex signal amplitude and noise
% variance.  Also returns a matched filtered output at symbol rate.
%
% INPUTS:
%   rx              - received signal at symbol rate
%   acqparam        - struct holding aquisition parameters
%   	.P              - position (indices) of pilot symbols
%     .p              - pilot symbols
%     .D              - position (indices) of data symbols
%     .M             - M for M-PSK, i.e size of constellation
%
% OUTPUTS:
%  ahat             - complex gain estimate
% Qhat             - value of objective function at estimator ahat (optional)
%
% (c) 2012 ITR-UniSA
% Authors: Robby McKilliam
% Created: 31 March 2012

    D = acqparams.D;
    P = acqparams.P;
    p = acqparams.p(:);
    M = acqparams.M;
    L = length(D) + length(P);
    
    %transpose stuff if required.
    y = rx(:);
    
    %some repeatedly used variables
    w = 2*pi/M;
    eta = exp(1i*w);
    nu = eta - 1.0;
    
    %setup structures
    phi = angle(y(D));
    u = w*round(phi/w);
    z = [phi - u, (1:length(D))'] ;
    g = y(D).*exp(-1i*u);
    %disp([ g(1), g(end) ]);
    
    %sort z and get permutation
    s = sortrows(z);
    sigma = @(k) s(mod(k-1, length(D)) + 1, 2);
    
    %now run the search loop
    Y = sum(y(P).*conj(p)) + sum(g);
    ahat = Y / L;
    Qhat = abs(Y)^2 / L;
    for( k = 1:(M+1)*length(D) )
        Y = Y + nu*g(sigma(k));
        g(sigma(k)) = eta*g(sigma(k));
        Q = abs(Y)^2 / L;
        if( Q > Qhat )
            Qhat = Q;
            ahat = Y / L;
        end
    end
    
end