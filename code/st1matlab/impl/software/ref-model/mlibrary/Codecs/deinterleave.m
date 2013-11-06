% Computes de-interleaved data according to interleaver matrix
%--------------------------------------------------------------------------
% deinterleave.m - This function is used to deinterleave incoming data
% according to the defined interleaver array
%
% y = deinterleave(x,permutation)
%
% Inputs:
%   x- interleaved data
%   permutation - interleaver array
% Outputs:
%   y - de-interleaved data
% Comments:
% 
%--------------------------------------------------------------------------
% Author: ITR
% Project: ASRP
%--------------------------------------------------------------------------
% Copyright 2012 : 
% Institute for Telecommunications Research
% University of South Australia
%--------------------------------------------------------------------------
% Contains ITR Background IP
%--------------------------------------------------------------------------

function y = deinterleave(x,permutation)

[m,n] = size(x);
switch length(permutation)
case m
  y(permutation,:) = x;
case n
  y(:,permutation) = x;
otherwise
  error('Length of permutation must be equal to one of the dimensions of x');
end
