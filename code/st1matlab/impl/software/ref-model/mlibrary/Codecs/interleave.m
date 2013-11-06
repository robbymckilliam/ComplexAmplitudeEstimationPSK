% Computes interleaved data according to interleaver matrix
%--------------------------------------------------------------------------
% interleave.m - This function is used to deinterleave incoming data
% according to the defined interleaver array
%
% y = interleave(x,permutation)
%
% Inputs:
%   x- data
%   permutation - interleaver array
% Outputs:
%   y - interleaved data
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

function y = interleave(x,permutation)

[m,n] = size(x);
switch length(permutation)
case m
  y = x(permutation,:);
case n
  y = x(:,permutation);
otherwise
  error('Length of permutaton must be equal to one of the dimensions of x');
end
