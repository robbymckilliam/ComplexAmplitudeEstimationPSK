function code = readldpccode(filename)
%readldpccode
%
%   reads LDPC struct from a file
%
%   code = readldpccode(filename)
%
%   code     ... ldpc code struct
%   filename ... file to read from
%
%   (c) 2006 Gottfried Lechner

fid = fopen(filename,'r');

N = fscanf(fid, '%i', 1);
M = fscanf(fid, '%i', 1);
E = fscanf(fid, '%i', 1);

code.vardegree   = fscanf(fid, '%i', N)';
code.chkdegree   = fscanf(fid, '%i', M)';
code.interleaver = fscanf(fid, '%i', E)';

code.interleaver = code.interleaver + 1;

fclose(fid);
