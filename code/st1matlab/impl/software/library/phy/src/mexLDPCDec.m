% mexLDPCDec   MATLAB wrapper for LDPC class
%    This mex-file enables the use of the LDPC class from MATLAB
%
%    Syntax:
%    The first argument is the ID of the LDPC code (starting from 0). This
%    allows to use many LDPC codes in parallel. The maximum number of
%    parallel codes is defined as MAXLDPC (default=100).
%
%    read the decoder:
%       mexLDPCDec(ID, 'readdecoder', filename);
%
%    encode:
%       parity = mexLDPCDec(ID, 'encodeLDGM', info);
%
%       Assumes that the code represents an LDGM code and computes parity bits.
%
%       parity = mexLDPCDec(ID, 'encodeRA', info);  
%
%       Assumes that the code represents an RA code and computes parity bits.
%
%    decode:
%       [Lapp, iterations, Isyn] = mexLDPCDec(ID, 'decode', Lch, maxit);
%
%       Performs decoding of the LLR-values Lch with a maximum number of maxit
%       iterations. The a-posteriori LLR-values Lapp and the number of
%       required iterations are returned.
%
%    code parameters:
%       [N, M] = mexLDPCDec(ID, 'dimensions');
%
%       Returns the size of the parity check matrix.
%
%    decoder configuration:
%       mexLDPCDec(ID, 'configure', parameter, value)
%
%       Possible parameters:
%         'clearmsg'
%           If true (default), the decoder clears the messages from check to variable nodes.
%           If false, the decoder uses the messages from the previous call. 
%         'method'
%           Defines the decoding method: 0 ... sum-product (default)
%                                        1 ... min-sum
%         'corrvec'
%           Sets the vector of correction factors for MSA decoding with post-processing
%
%    status information:
%       b   = mexLDPCDec(ID, 'decloaded')
%
%   AUTHOR: Gottfried Lechner (gottfried.lechner@unisa.edu.au)
%   Copyright 2010-2012 Gottfried Lechner ITR-UniSA
%