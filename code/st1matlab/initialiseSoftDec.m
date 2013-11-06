% Initialises SoftDec structure
%---------------------------------------------------------------------
% SoftDec = initialiseSoftDec(PacketLength, RxLength, SatChanParams)
%
% Inputs:
%   PacketLength    number of symbols
%   RxLength        number of samples
%   SatChanParams   channel parameters (optional)
%
% Outputs:
%   SoftDec         empty SoftDec structure
% 
%---------------------------------------------------------------------
% Author: Gottfried Lechner
% Project: ASRP
%---------------------------------------------------------------------
% Copyright 2013 
% Institute for Telecommunications Research
% University of South Australia
%---------------------------------------------------------------------

function SoftDec = initialiseSoftDec(PacketLength, RxLength, SatChanParams)

SoftDec.Fo        = [];
SoftDec.Fd        = [];
SoftDec.FdDir     = 1;
SoftDec.Td        = [];
SoftDec.PhOff     = [];
SoftDec.SigPowEst = 0;
SoftDec.NoiseVar  = [];
SoftDec.Acquired  = 0;
SoftDec.success   = 0;
SoftDec.finish    = 0;
SoftDec.SoftSyms  = zeros(1,PacketLength);
SoftDec.RxMsgDec  = [];
SoftDec.RxSyms    = [];
SoftDec.SampHat   = zeros(1,RxLength);
SoftDec.SampHatp  = [];
SoftDec.TTL       = 0;
SoftDec.DecStat   = 0;
SoftDec.DecStatP  = 0;
SoftDec.Fs        = 0;

if nargin>2
  SoftDec.Acquired  = 1;
  SoftDec.Fo        = SatChanParams.Fo;
  SoftDec.Td        = SatChanParams.Td;
  SoftDec.PhOff     = SatChanParams.PhOff;
  SoftDec.Fd        = SatChanParams.Fd;
  SoftDec.FdDir     = SatChanParams.FdDir;
  SoftDec.SigPowEst = SatChanParams.Esam*10^(SatChanParams.Power/10);   
end
  
end
