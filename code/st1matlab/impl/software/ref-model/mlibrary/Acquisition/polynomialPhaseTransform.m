function [phiShifted,freqShifted,frateShifted] = polynomialPhaseTransform(phi,freq,frate,tau)

% This function applies a polynomial phase transform to change the time
% reference point (t=0). The user-provided coefficients describe the
% second-order polynomial phase signal
%
%       phi + 2*pi*freq*t + 2*pi*(frate/2)*t^2,
%
% which are used to calculate the coefficients of the polynomial
%
%       phiShifted + 2*pi*freqShifted*(t-tau) + 2*pi*(frateShifted/2)*(t-tau)^2.
%
% Relative to the original phase polynomial, the time origin of the new
% polynomial is shifted by tau.
%
% For details, see Andre Pollok's ASRP workbook 007, page 39.
%
% (c) 2012 ITR-UniSA
% Authors: Andre Pollok
% Created: 29 July 2012

frateShifted = frate;
freqShifted  = freq + 2*(frate/2)*tau;
phiShifted   = phi  + 2*pi*freq*tau + 2*pi*(frate/2)*tau^2;