% Reads parity-check matrix from PEG file format in compressed matrix format
%---------------------------------------------------------------------
% readsparsepegmatrix.m - reads parity-check matrix from PEG file format
% in compressed matrix format
% 
%   c = readsparsepegmatrix(filename)
%
% Inputs:
%   filename    - Filename of the PEG file
% Outputs:
%    c          - LDPC code in compressed matrix format
% Comments:
%      
%---------------------------------------------------------------------
% Author:   G. Lechner, ITR
% Project : ASRP
%---------------------------------------------------------------------
% Copyright 2011 : 
% Institute for Telecommunications Research
% University of South Australia
%---------------------------------------------------------------------
%
%---------------------------------------------------------------------

function c = readsparsepegmatrix(filename)


fid = fopen(filename, 'r');

N   = fscanf(fid, '%i',1);
M   = fscanf(fid, '%i',1);
dc  = fscanf(fid, '%i',1);

sparsepeg = zeros(M,dc);
for m=1:M
  for d=1:dc
    sparsepeg(m,d) = fscanf(fid, '%i',1);
  end
end
fclose(fid);

[dv,x]= hist(nonzeros(sparsepeg(:)),1:N);
dvmax = max(dv);
c     = zeros(N,dvmax);

for n=1:N
  idx = find(sum((sparsepeg==n),2))';
  idx = [idx zeros(1,dvmax-size(idx,2))];

  c(n,:) = idx;
end
