% Calls encoder function
%---------------------------------------------------------------------
% encoder.m - Encodes data using the LDPC Codec definition
%
%   [TxEnc] = encoder(TxUncData)
%
% TBD
%
% Inputs:
%   TxUncData   -   Binary uncoded data
%   UserParams  -   User Codec Information
% Outputs:
%   TxEnc       -   Coded Data
% Comments:
% 
%---------------------------------------------------------------------
% Author:  G. Lechner, ITR
% Project: ASRP
%---------------------------------------------------------------------
% Copyright 2011 : 
% Institute for Telecommunications Research
% University of South Australia
%---------------------------------------------------------------------
% 
%---------------------------------------------------------------------

function  [TxEnc] =encoder(TxUncData,UserParams)

      TxEnc=rem(TxUncData*UserParams.Hu',2);
      
      if strcmp(UserParams.ModScheme, 'qpsk')
        TxEnc = [TxUncData rem(cumsum(TxEnc),2)];
      end

end