% Writes LDPC struct to a file
%---------------------------------------------------------------------
% writeldpccode.m - writes LDPC struct to a file
% 
% writeldpccode(code, filename)
%
% Inputs:
%   code     - LDPC code struct
%   filename - File to write the code
% Outputs:
%   None
% Comments:
%      
%---------------------------------------------------------------------
% Author:   ITR
% Project : ASRP
%---------------------------------------------------------------------
% Copyright 2006 : 
% Institute for Telecommunications Research
% University of South Australia
%---------------------------------------------------------------------
% Contains ITR Background IP
%---------------------------------------------------------------------
function writeldpccode(code, filename)


N   = length(code.vardegree);
M   = length(code.chkdegree);
E   = length(code.interleaver);

fid = fopen(filename,'w');

fprintf(fid, '%i %i %i\r\n', N, M, E);
fprintf(fid, '%i ', code.vardegree);
fprintf(fid, '\n');
fprintf(fid, '%i ', code.chkdegree);
fprintf(fid, '\n');
fprintf(fid, '%i ', code.interleaver-1);
fprintf(fid, '\n');

fclose(fid);
