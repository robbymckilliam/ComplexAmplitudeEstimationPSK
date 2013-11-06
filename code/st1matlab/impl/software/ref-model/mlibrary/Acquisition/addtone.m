% Adds a Tone to main signal to assist ST1/ST2 acquisition  
%--------------------------------------------------------------------------
% addtone - Creates and adds a sine tone to the main Tx samples
%
% [ OutTxSamples ] = addtone( InTxSamples,UserParams,ModemParams,STAcqParams )
%
% Inputs:
%   InTxSamples - Main Tx samples
%   UserParams  - Structure containg user parameters 
%   ModemParams - Structure containing simulation parameters
%   STAcqParams - Structure used by ST1/ST2 to add tone to tx signal
% Outputs:  
%   OutTxSamples - Signal containing main tx samples and added sine tone
% Comments:
% 
%--------------------------------------------------------------------------
% Author:  
% Project: ASRP
%--------------------------------------------------------------------------
% Copyright 2012 : 
% Institute for Telecommunications Research
% University of South Australia
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function OutTxSamples = addtone(InTxSamples,~,ModemParams,STAcqParams,scaling)

% if no scaling is provided then don't scale
if nargin<5
  scaling = 1;
end

%-------------------------------------------------------------------------%
%   Create sinusoidal tone
%-------------------------------------------------------------------------%
% Lsin is the number of symbols spanned by the tone
% Hence, the length of the tone is Lsin-1 symbol durations plus leading and
% trailing filter tails.
ToneLenSamples  = (STAcqParams.Lsin-1)*ModemParams.P + 2*ModemParams.RrcNsym*ModemParams.P+1;

% set up indices for tone and set zero phase of tone to coincide with the
% first symbol (instead of the first sample) to match hardware
% implementation (Andre, Gottfried, John)
SampleTime  = (0:(ToneLenSamples-1))/ModemParams.Fs;
SinDelay    = ((length(ModemParams.TxRrcCoefs)-1)/2+1)/ModemParams.Fs;

% create the tone with constant amplitute of one
cexp = zeros(1,ToneLenSamples);
cexp(1:(ToneLenSamples)) = STAcqParams.amplsin*...
  exp(2i*pi*STAcqParams.fsin*(SampleTime-SinDelay));

%-------------------------------------------------------------------------%
%   Create ramp
%-------------------------------------------------------------------------%
% sample indices for actual ramp
RampUpIndex = 1:(STAcqParams.ToneRampSym*ModemParams.P);

% weights for actual ramp
ramp        = (sin(pi/2 *...
  ((RampUpIndex-1)/ModemParams.Fs) / ...
  (STAcqParams.ToneRampSym/ModemParams.SymbolRate))).^2;

% prepend half a symbol-period of zeros
ramp        = [zeros(1,ModemParams.P/2) ramp];

% set up window for the duration of the tone
window  = ones(size(cexp));
% copy ramp to start of window
window(1:length(ramp))              = ramp;
% copy mirrored ramp to end of window (assuming symmetric up/down ramps)
window(((end-length(ramp)+1):end))  = fliplr(ramp);

%-------------------------------------------------------------------------%
%   Add tone to Tx signal (apply window and scaling)
%-------------------------------------------------------------------------%           
OutTxSamples = InTxSamples;
OutTxSamples(1:ToneLenSamples) = InTxSamples(1:ToneLenSamples) + scaling*cexp.*window;

end

