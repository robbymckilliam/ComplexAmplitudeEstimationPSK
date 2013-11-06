% Reads LDPC structure from a file
%---------------------------------------------------------------------
% readldpccode.m - Reads LDPC structure from a file
% 
% code = readldpccode(filename)
%
% Inputs:
%   filename    - File to read from
% Outputs:
%    code       - LDPC code structure
% Comments:
%      
%---------------------------------------------------------------------
% Author:   G.Lechner,ITR
% Project : ASRP
%---------------------------------------------------------------------
% Copyright 2006 : 
% Institute for Telecommunications Research
% University of South Australia
%---------------------------------------------------------------------
% Contains ITR Background IP
%---------------------------------------------------------------------
function code = readldpccode(filename)

fid = fopen(filename,'r');

N = fscanf(fid, '%i', 1);
M = fscanf(fid, '%i', 1);
E = fscanf(fid, '%i', 1);

code.vardegree   = fscanf(fid, '%i', N)';
code.chkdegree   = fscanf(fid, '%i', M)';
code.interleaver = fscanf(fid, '%i', E)';

code.interleaver = code.interleaver + 1;

fclose(fid);
