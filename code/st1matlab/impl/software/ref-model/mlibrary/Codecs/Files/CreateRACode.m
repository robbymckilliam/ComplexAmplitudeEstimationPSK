% Script to construct a parity-check matrix for RA coded modulation.
%---------------------------------------------------------------------
% CreateRACode.m - Example script to construct a parity-check matrix
%   for RA coded modulation.
% 
% 
% CreateRACode
%
% Inputs:
%   None
% Outputs:
%   None
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

datalen= 440; % number of input bits into R1/2 codec
vardegrees  = 'RA1.vdeg';
chkdegrees  = 'RA1.cdeg';
peg         = '../PEG/peg';
name        = ['RA1N' num2str(datalen)];
N           = datalen*2;
M           = datalen;

% build commandline for PEG and run to create LDGM part
cmdline     = sprintf('%s -numM %i -numN %i -codeName %s.peg -degFileName %s -degFileNameChk %s -sglConcent 0', peg, M, N-M, name, vardegrees, chkdegrees);
system(cmdline);

% read PEG file
c   = readsparsepegmatrix(sprintf('%s.peg', name));

% convert to sparse matrix
Hu  = compressed2matrix(c);

% interleave the matrix to remove structure of the PEG algorithm
Hu  = Hu(randperm(M),:);
Hu  = Hu(:,randperm(N-M));

% add accumulator matrix
A   = speye(M);
A   = A + [zeros(1,M) ; A(1:(M-1),:)];
H   = [Hu A];

% convert to LDPC code structure
code = matrix2ldpccode(H);

% write LDPC code to file
writeldpccode(code, sprintf('%s.dec', name));

% cleanup
delete(sprintf('%s.peg.lhg', name));
delete(sprintf('%s.peg', name));
