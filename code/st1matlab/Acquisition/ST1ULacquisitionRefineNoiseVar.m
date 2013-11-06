% Refine noise variance estimate.
%-------------------------------------------------------------
% ST1ULacquisitionRefineNoiseVar.m 
%
% [SoftDec] = ST1ULacquisitionRefineNoiseVar(Rx, SoftDec, ModemParams, STAcqParams)
% 
% This function should be called after ST1ULacquisitionRefinePeriodoWrapper
% and refines the noise variance estimate in code-aided fashion, i.e. the
% soft symbols provided by the decoder are treated as pilot symbols.
% 
% Inputs:
%   Rx              - Vector containing samples of received signal
%   SoftDec         - ST1UL SoftDec struct
%   ModemParams     - ST1UL ModemParams struct
%   AcqParams       - Struct holding aquisition parameters
%       .P              - position (indices) of pilot symbols
%       .p              - pilot symbols
%       .D              - position (indices) of data symbols
%       .M              - M for M-PSK, i.e. size of constellation
%       ToneInSymbols   - Vector containing the tone after matched 
%                         filtering and downsampling to symbol rate.
%
% Outputs:
%   SoftDec         - ST1UL SoftDec struct with refined noise variance
%                     estimate
%
% Comments:
% 
%-------------------------------------------------------------
% Author: A. Pollok, Institute for Telecommunications Research
% Project: ASRP
%-------------------------------------------------------------
% Copyright 2013
% Institute for Telecommunications Research
% University of South Australia
%-------------------------------------------------------------

function [SoftDec] = ST1ULacquisitionRefineNoiseVar(Rx, SoftDec, ModemParams, STAcqParams)

% get matched filter output
ModemParams.ChanType    = 'chancomp';
DataComp                = satchan(Rx, SoftDec, ModemParams, SoftDec, [],  STAcqParams);
MFoutput                = DataComp/sqrt(SoftDec.SigPowEst); %#ok<*AGROW>
MFoutput                = MFoutput(1:ModemParams.CodedSymLenR);

% tone cancellation
if STAcqParams.CancelTone
    MFoutput = MFoutput - STAcqParams.ToneInSymbols;
end

% re-estimate Es/No, treating decoder soft symbols as pilots
STAcqParams.P = [STAcqParams.P STAcqParams.D];
STAcqParams.p = SoftDec.SoftSyms;
STAcqParams.D = [];
EsNohat = EsNoest(MFoutput, STAcqParams);
SoftDec.NoiseVar = SoftDec.SigPowEst/EsNohat;

end