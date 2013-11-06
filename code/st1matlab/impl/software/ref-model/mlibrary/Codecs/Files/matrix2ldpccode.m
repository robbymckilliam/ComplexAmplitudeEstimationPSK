% Converts sparse parity-check matrix to LDPC code struct
%---------------------------------------------------------------------
% matrix2ldpccode.m -  Converts sparse parity-check matrix to LDPC 
% code struct
%
%  code = matrix2ldpccode(H)
% 
%
% Inputs:
%    H          - sparse parity-check matrix
% Outputs:
%    code       - LDPC code struct
% Comments:
%      
%---------------------------------------------------------------------
% Author:   G. Lechner, ITR
% Project : ASRP
%---------------------------------------------------------------------
% Copyright 2010 : 
% Institute for Telecommunications Research
% University of South Australia
%---------------------------------------------------------------------
% Contains ITR Background IP
%---------------------------------------------------------------------
function code = matrix2ldpccode(H)

c = zeros(size(H,2), max(sum(H,1)));
for j=1:size(H,2)
   x = find(H(:,j));
   c(j,1:length(x)) = x; 
end

N = size(c,1);
M = max(max(c));

% vardegree
code.vardegree = sum((c>0),2)';

% chkdegree
x = zeros(1,M);
for i=1:N
    idx = nonzeros(c(i,:));
    x(idx) = x(idx) + 1;
end
code.chkdegree = x;

% interleaver
y = sort(c');
y = nonzeros(y(:));
x = zeros(1,sum(code.vardegree));
j = 1;
for i=1:M
    idx = find(y==i);
    s   = size(idx,1);
    x(j:j+s-1) = idx;
    j = j+s;
end
code.interleaver = x;

end
