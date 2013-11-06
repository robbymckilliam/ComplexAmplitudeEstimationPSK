% Generates a B-ary discretee memoryless source
%---------------------------------------------------------------------
% DMS.m - B-ary discrete memoryless source, returns L independently 
% chosen uniformly distributed symbols in the range [0:B-1]
%
%   b = DMS(L,B)
%
% TBD
%
% Inputs:
%  L - length of data
%  B - order
% Outputs:
%   b - discrete memoryless B-ary data
% Comments:
%   Care should be taken with seeding rand() so that repeated
%   simulations do not get the same data
%---------------------------------------------------------------------
% Author:  ITR
% Project: ASRP
%---------------------------------------------------------------------
% Copyright 2011 : 
% Institute for Telecommunications Research
% University of South Australia
%---------------------------------------------------------------------
% Contains ITR Background IP
%---------------------------------------------------------------------
function b = DMS(L,B)
    b = floor(B*(rand(1,L)));
end