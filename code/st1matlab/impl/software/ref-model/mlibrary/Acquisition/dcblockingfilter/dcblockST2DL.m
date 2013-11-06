% DC blocking filter
%-------------------------------------------------------------
% dcblock.m 
%
% [RxSampFilt] = dcblockST2DL(RxSamp)
%
% This digital filter blocks the DC component in the received signal. It is
% primarily intended to remove the complex exponential at DC that is
% transmitted in the ST2 downlink waveform for frequency estimation
% purposes. Once the frequency has been estimated, the discrete tone should
% be removed in order to avoid bit error rate degradations. The
% implementation uses an order-one IIR filter.
% 
% Inputs:
%   RxSamp      - Vector containing samples of received signal
%
% Outputs:
%   RxSampFilt  - Vector containing filtered received samples
%
% Comments: For details see [J. de Freitas, "The DC Blocking Filter", Jan.
%           2007, available at
%           http://www.mathworks.com.au/matlabcentral/fileexchange/13792-the-dc-blocking-filter] 
% 
%-------------------------------------------------------------
% Author: A. Pollok, Institute for Telecommunications Research
% Project: ASRP
%-------------------------------------------------------------
% Copyright 2012 
% Institute for Telecommunications Research
% University of South Australia
%-------------------------------------------------------------

function [RxSampFilt] = dcblockST2DL(RxSamp)

% ST1 UL: coefficents of order-one IIR filter

p = 0.998009283733233; % coeff for cut-on frequency at 12Hz
    % fc = 12 / (OS/(2*Ts)); % normalised cut-on frequency
    % p = dcblock_MatlabCentral(fc);
b = [1 -1];
a = [1 -p];

% filter frequency corrected samples
RxSampFilt = filter(b,a,RxSamp);

end