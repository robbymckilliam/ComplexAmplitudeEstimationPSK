% Joint pilot and data based Es/No estimator.
%-------------------------------------------------------------
% EsNoest.m 
%
% [EsNohat] = EsNoest(MFoutput,AcqParams)
%
% This function etimates the Es/No from the complex matched filter (MF)
% output. First, hard decisions for the data symbols are obtained from the
% MF output. After inserting the known pilot symbols, the symbol stream is
% used to strip the modulation of the MF output. The Es/No is estimated
% from the resulting signal.
% 
% Inputs:
%   MFoutput    - Vector containing complex matched filter output at one 
%                 sample per symbol.
%   AcqParams   - Struct holding aquisition parameters
%       .P              - position (indices) of pilot symbols
%       .p              - pilot symbols
%       .D              - position (indices) of data symbols
%       .M              - M for M-PSK, i.e. size of constellation
%
% Outputs:
%   EsNohat     - Estimate of the Es/No.
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

function [EsNohat] = EsNoest_coder(MFoutput,M,D,P,p)

symbhat = zeros(1,length(P)+length(D))+1i*zeros(1,length(P)+length(D));
% make hard decisions for PSK data symbols
symbsectorhat = mod(round(angle(MFoutput(D))/(2*pi/M)), M);
symbhat(D) = exp(2i*pi/M*symbsectorhat);
% insert pilot symbols
symbhat(P) = p;

% Bill's blind Es/No estimator
MFoutDemod = MFoutput .* conj(symbhat); % strip modulation off
EsNohat = abs(mean(MFoutDemod) / std(MFoutDemod))^2;

end