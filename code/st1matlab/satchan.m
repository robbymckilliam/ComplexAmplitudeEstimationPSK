% Adds or compensates for channel offset impairments in the tx signal
%--------------------------------------------------------------------------
% SatChan.m - adds or removes channel offsets (estimates) from incoming
% signal. This function runs in two modes, either to add offsets to
% to a multi-user array (each user gets their own offsets), in compensation
% mode (each user has offset estimates removed from signal)
%
% Inputs:
%   Syms - transmitted symbol sequence (can be multi-user array)
%   InChanParams - Channel Offset structure
%   ModemParams - Main modem parameters structure
%   UserParams  - User parameters structure
%   SatChanMode - optional, if defined in function call, AWGN is added to output signal
% Outputs:
%   RxSyms - Symbol sequence with channel offsets added/compensated
%   OutChanParams - Updated channel offsets (in case of Doppler)
%   OutUserParams - Debug only (to store data snaphots)
%
%
%--------------------------------------------------------------------------
% Author:  ITR
% Project: ASRP
%--------------------------------------------------------------------------
% Copyright 2011 :
% Institute for Telecommunications Research
% University of South Australia
%--------------------------------------------------------------------------
% Contains ITR Background IP
%--------------------------------------------------------------------------
function [RxSyms, OutChanParams,OutUserParams] = satchan(...
  SymsIn, InChanParams, ModemParams, UserParams, SatChanMode, STAcqParams) %#ok<INUSD>


if nargin<6
  STAcqParams = [];
end

if nargin<5
  SatChanMode = [];
end

OutChanParams = InChanParams;
OutUserParams = UserParams;

interpolation = 'linear';

%-------------------------------------------------------------------------%
%   Check Size of Incoming Data (can be multi-dimensional or single user)
%-------------------------------------------------------------------------%
[K,N]   = size(SymsIn);
%-------------------------------------------------------------------------%
%   Assign Array for full length of Slot (inc Guard Time)
%-------------------------------------------------------------------------%
ChanOut = zeros(1,ModemParams.SlotLenSam);

%-------------------------------------------------------------------------%
%  Add Channel Offsets to the Signal ('chanadd' mode)
%-------------------------------------------------------------------------%
if strcmp(ModemParams.ChanType,'chanadd') % add offsets mode
  for k=1:K
    SymsPacket  = zeros(1,length(ChanOut)); %user slot storage
    
    %-----------------------------------------------------------------%
    %   Add Integer Timing Offset (in number of samples)
    %-----------------------------------------------------------------%
    SymsPacket(1,round(InChanParams(k).Td)+1:round(InChanParams(k).Td)+N) = SymsIn(k,:);
    
    %-----------------------------------------------------------------%
    %   Add Fractional Timing Offset
    %-----------------------------------------------------------------%
    x           = 1:length(SymsPacket);
    xi          = x-(InChanParams(k).Td - round(InChanParams(k).Td)); %apply sampling either early(-ve) or late(+ve)
    SymsPacket  = interp1(x,SymsPacket,xi, interpolation,0);
    
    %-----------------------------------------------------------------%
    %   Add phase (note factor of 2 in Doppler (f(t)=f+2*Dop*t))
    %-----------------------------------------------------------------%
%     time        = 0:(ModemParams.SlotLenSam-1);
    time        = 0:(length(SymsPacket)-1);
    Theta       = InChanParams(k).PhOff +...
                  2*pi*InChanParams(k).Fo*time +...
                  2*pi/2*InChanParams(k).FdDir*InChanParams(k).Fd/ModemParams.Fs*time.^2;
    SymsPacket  = SymsPacket.*exp(1j*Theta);
    
    %-----------------------------------------------------------------%
    %   Handle Signal Scaling
    %-----------------------------------------------------------------%
    if ~isempty(SatChanMode)
      %scale incoming signal before adding noise
      SymsPacket = SymsPacket.*sqrt(InChanParams(k).Esam);
      % assign power to individual packets
      SymsPacket = SymsPacket.*10^(InChanParams(k).Power/20);
    else
      % scale with signal amplitude estimate
      % Note - in perfect chan mode, this is the know amp+power set
      % in snrest.m
      SymsPacket = SymsPacket.*sqrt(InChanParams(k).SigPowEst);
    end
    
    if ~isempty(SatChanMode) %DEBUG storage
      OutUserParams(k).TxSymsChan = SymsPacket;
    end
    %-------- Create combined users signal----------------------------%
    if length(ChanOut)<length(SymsPacket)
      ChanOut     = [ChanOut zeros(1,length(SymsPacket)-length(ChanOut))];
    else
      SymsPacket  = [SymsPacket zeros(1,length(ChanOut)-length(SymsPacket))];
    end
    ChanOut = ChanOut + SymsPacket; %accumulate user signals    
  end
  %---------------------------------------------------------------------%
  %   Add AWGN to Combined Users Signal
  %---------------------------------------------------------------------%
  if nargin>4
    % Add AWGN to scaled signal, where N0=1
    RxSyms = ChanOut + sqrt(1/2)*(randn(1,length(ChanOut)) + 1j*randn(1,length(ChanOut)));
  else
    % channel rebuild, no noise added, simply output scaled packet
    RxSyms = ChanOut;
  end
  RxSyms = RxSyms(1:ModemParams.SlotLenSam);
end

%-------------------------------------------------------------------------%
%  Remove Channel Offsets from the Signal ('chancomp' mode)
%-------------------------------------------------------------------------%
if strcmp(ModemParams.ChanType,'chancomp') % add offsets mode
  % Remove phase
  time    = 0:(length(SymsIn)-1);
  Theta   = InChanParams.PhOff +...
            2*pi*time*InChanParams.Fo +...
            2*pi/2*InChanParams.FdDir*InChanParams.Fd/ModemParams.Fs*time.^2;
  RxSyms  = SymsIn.*exp(-1j*Theta);
  
  %---------------------------------------------------------------------%
  %     LPF signal prior to fractional interpolation
  %---------------------------------------------------------------------%
  lpf_coefs_len = 21; % length of LPF for timing interpolation
  lpf_coefs     = fir1(lpf_coefs_len-1,0.5);      % create LPF filter
  RxSyms        = upfirdn(RxSyms,lpf_coefs,1,1);
  RxSyms(1:(lpf_coefs_len-1)/2)         = [];     % chop off filter tails
  RxSyms(end-(lpf_coefs_len-1)/2+1:end) = [];     % chop off filter tails
  
  %---------------------------------------------------------------------%
  %     Remove Fractional Timing Offset
  %---------------------------------------------------------------------%
  x       = 1:length(RxSyms);
  xi      = x+(InChanParams.Td - round(InChanParams.Td));
  RxSyms  = interp1(x,RxSyms,xi, interpolation,0);
  
  %---------------------------------------------------------------------%
  %     Remove Integer Timing Offset
  %---------------------------------------------------------------------%
  Nb      = N-ModemParams.TguardSam;
  RxSyms  = RxSyms(round(InChanParams.Td)+1:round(InChanParams.Td)+Nb);
  
  % filter tone
  % DEBUG: need to work out the influence of cancelling the tone on the SNR
  % estimator. Dodgy fix at the moment.
  if ~isempty(STAcqParams)
    if STAcqParams.FilterTone
      ts      = (0:(length(RxSyms)-1))/ModemParams.Fs;
      Ebefore = mean(abs(RxSyms.^2));
      RxSyms  = RxSyms .* exp(-2i*pi*STAcqParams.fsin*ts);
      RxSyms  = dcblockST1UL(RxSyms);
      RxSyms  = RxSyms .* exp(+2i*pi*STAcqParams.fsin*ts);
      Eafter  = mean(abs(RxSyms.^2));
      RxSyms  = RxSyms * sqrt(Ebefore/Eafter);
    end
  end
  
  %-------------------------------------------------------------%
  %     RRC Rx Pulse Shaping filter
  %-------------------------------------------------------------%
  Lx      = length(RxSyms);
  RxSyms  = upfirdn(RxSyms,ModemParams.RxRrcCoefs,1,1);
  Ly      = length(RxSyms);
  Lt      = Ly-Lx;
  Lht     = length(ModemParams.TxRrcCoefs);
  Lhr     = length(ModemParams.RxRrcCoefs);
  RxSyms(1:(Lhr-1)/2+(Lht-1)/2)    = [];%chop front Tx+Rx tails
  %RxSyms(end-(Lt-(Lhr-1)/2)+1:end) = [];%chop end Rx tail
  RxSyms(end-((Lhr-1)/2+(Lht-1)/2)+1:end) = [];%chop end Rx tail
  
  %-------------------------------------------------------------%
  %     Downsample
  %-------------------------------------------------------------%
  RxSyms  = RxSyms(1:ModemParams.P:length(RxSyms));    
end
  
end


