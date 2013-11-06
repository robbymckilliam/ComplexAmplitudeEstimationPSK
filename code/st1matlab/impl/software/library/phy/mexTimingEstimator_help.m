% mexTimingEstimator   MATLAB wrapper for TimingEstimator class
%
%   This mex-file enables the use of the TimingEstimator class from MATLAB.
%
%   Construct a ST1 timing estimator 
%       mexTimingEstimator(name,P,D,pilots,T,Ts,p,q,taumin,taumax,rolloff,duration,tol);
%
%   Construct a ST2 timing estimator 
%       mexTimingEstimator(name,P,D,pilots,T,Ts,taumin,taumax,rolloff,duration);
%   
%   Run timing estimator
%       tauhat = mexTimingEstimator(name,r);
%
%    Input parameters:
%       name      - name of this estimator, a string. 
%       P           - position (indices) of pilot symbols
%       D           - position (indices) of data symbols       
%       pilots      - complex pilot symbols
%       T           - symbol period
%       Ts         - sample period
%       p, q       - relatively prime integers such that Ts=p*T/q.
%                     For example for 4 times oversampling p = 1
%                     and q = 4.
%       taumin      - minimum time offset
%       taumax      - maximum time offset
%       rolloff     - rolloff for the root raised cosine transmit pulse 
%       duration - duration of the transmit pulse
%       tol         - [OPTIONAL] tolerance for time estimate. If omitted,
%                     a default value of 1e-6 is used.
%       r            - received complex signal
%
%   Output parameters:
%       tauhat      - timing offset estimate
%
%   (C) 2012 Andre Pollok, ITR/UniSA
%