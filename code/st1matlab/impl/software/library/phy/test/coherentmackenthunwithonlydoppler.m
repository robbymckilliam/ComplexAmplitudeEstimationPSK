function [ahat,fhat] = coherentmackenthunwithonlydoppler(rx, acqparams)
% Runs the phase estimator a number of times at different Doppler
% and Doppler rate.  Returns the estimated gain, amplitude, Doppler
% and Doppler rate.
    
% INPUTS:
%   rx              - received signal at symbol rate
%   acqparam        - struct holding aquisition parameters
%       .P              - position (indices) of pilot symbols
%       .D              - position (indices) of data symbols
%     .p              - pilot symbols
%     .M             - M for M-PSK, i.e size of constellation
%     .T              - symbol period
%     .f_max      - maximum Doppler
%     .f_min       - minimum Doppler
%
% OUTPUTS:
%  ahat             - complex gain estimate
%
% (c) 2012 ITR-UniSA
% Authors: Robby McKilliam, Andre Pollok
% Created: 21 May 2012, Robby McKilliam
% Modified: 25 May 2012, Andre Pollok

    M = acqparams.M;
    D = acqparams.D;
    P = acqparams.P;
    p = acqparams.p(:);
    L = length(D) + length(P);
    T = acqparams.T; % symbol period
    max_index = max([P,D]);
    inddif = max_index - min([P,D]);

    f_min = acqparams.f_min;
    f_max = acqparams.f_max;
    
    SEARCH_OVERSAMPLE = 2.0; %adjusts how fine the search grid is
                           %(smaller is faster, larger is more accurate)
    
    f_trials = ceil(SEARCH_OVERSAMPLE*M*(T*inddif)*(f_max - f_min));
    fstep = (f_max - f_min)/f_trials;
    if( isnan(fstep) ) 
        fstep = 1;
    end
    
    %symbol times vector
    ts = T*[P,D];
     
    %some repeatedly used variables
    w = 2*pi/M;
    eta = exp(1i*w);
    nu = eta - 1.0;
    
    bestQ = 0;
    fhatc = f_min;
    ahatc = 0;
    for f = f_min:fstep:f_max
            
            %compute the Doppler rate shifted recieved signal.
            y = zeros(1,max_index);
            y([P,D]) = rx([P,D]) .* exp(-2i*pi*(f*ts));
            
            %run the least squares estimator on the shifted signal
            [a, Q] = coherentmackenthun(y, acqparams);
            
            if( Q > bestQ )
                bestQ = Q;
                ahatc = a;
                fhatc = f;
            end
            
       
    end
    
    %now setup an optimiser.  Will use matlab's fminbnd, but there
    %are stubs here for Newton's method
    
    %first construct the vector s with the hard decisions and pilots
    s = zeros(1, max_index);
    s(P) = p; %set the pilots
    u = round(M/(2*pi)*(angle(rx(D)) - angle(ahatc) - 2*pi*(fhatc*T*D)));
    s(D) = exp(2i*pi/M * u ); %set the hard decisions
    
    %get the recieved signal at the pilot positions
    y = rx([P,D]);
      
    %% refine estimates via Newton's method
    % For details, see Andre Pollok's ASRP workbook, page 34-35.
    
    alpha = 0; % alpha=0: max log(periodogram), alpha=1: max standard periodogram
    newtonstep = Inf(2,1);
    EPSILON = 1e-10;
    MAX_ITER = 10;
    numIter = 0;
    fhat = fhatc;
    r = y .* conj(s([P,D]));
    while ( any(abs(newtonstep) > EPSILON) && (numIter < MAX_ITER))
        % evaluate periodogram I and its first two derivatives Idd and Idd
        X = r .* exp(-2i*pi * (fhat*ts));
        Y = 1/L*sum(X);
        Ydf     =  -2i*pi/L * sum(ts   .*X);
        Ydf2    = -4*pi^2/L * sum(ts.^2.*X);
        I       = abs(Y)^2;
        Idf     = 2*real( Ydf * conj(Y) );
        Idf2    = 2*real( Ydf2 * conj(Y)) + 2*abs(Ydf)^2;
        
        G = Idf^2;
        H = Idf2;
        
        % Newton iterate with monotonic function kappa(I) = I^alpha.
        % For alpha=1, this is the standard Newton iterate inv(H)*Idf.
        newtonstep = Idf / ((alpha-1)/I*G + H);
        fhat   = fhat   - newtonstep(1);

        numIter = numIter + 1;
    end
    
        X = r .* exp(-2i*pi * (fhat*ts));
     Y = 1/L*sum(X);
    ahat = Y;
    %ahat = ahatc;

end