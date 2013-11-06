% Converts parity-check matrix from compressed matrix to sparse matlab matrix
%---------------------------------------------------------------------
% compressed2matrix.m - converts parity-check matrix from compressed 
% matrix to sparse matlab matrix
% 
% H = compressed2matrix(c)
%
% Inputs:
%    c           - compressed matrix format
% Outputs:
%    H           - sparse parity-check matrix
% Comments:
%      
%---------------------------------------------------------------------
% Author:   ITR
% Project : ASRP
%---------------------------------------------------------------------
% Copyright 2011 : 
% Institute for Telecommunications Research
% University of South Australia
%---------------------------------------------------------------------
% 
%---------------------------------------------------------------------
function H = compressed2matrix(c)

rows = max(max(c));
cols = size(c,1);

H = spalloc(rows,cols,length(nonzeros(c)));

for i=1:cols
   H(nonzeros(c(i,:)),i) = 1;
end
