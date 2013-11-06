function [ModemParams, STAcqParams] = overwriteDefaultParams(Modem, STAcq)

%% Modem Parameters
ModemParams.Codec         = 'ST1UP_MCS0_RA1N232.dec';
ModemParams.InfoLen       = 232;
ModemParams.CodewordLen   = 464;
ModemParams.S             = st1interleaver();
ModemParams.K             = [];
ModemParams.Hu            = [];

ModemParams.SymbolRate    = 1024;
ModemParams.P             = 32;
ModemParams.Fs            = [];

ModemParams.ModScheme     = 'qpsk';
ModemParams.AddTone       = 1;
ModemParams.ModulationMap = [-1-1j -1+1j 1-1j 1+1j]/sqrt(2);
ModemParams.RrcAlpha      = 0.5;
ModemParams.RrcNsym       = 3.5;  % half length of the pulse shape filter
ModemParams.L             = 0;
ModemParams.AISFormat     = 0;

ModemParams.N_pilots      = 1;
ModemParams.pilot.pos     = 1;
ModemParams.pilot.bits    = [0 0 1 1 0 0 1 1 0 0 0 0 1 1 0 1 1 0 0 1 1 0 0 1 0 1 1 0];
ModemParams.pilotindex    = [];

ModemParams.SoftDecoders  = 5;
ModemParams.MUDIter       = 50;
ModemParams.TTL           = 50;
ModemParams.maxit         = 50;
ModemParams.RefineIter    = 2;
ModemParams.PerfectChan   = 0;

ModemParams.TguardSam     = 0;

% update default structure with user supplied values
if nargin>=1
  if ~isempty(Modem)
    ModemParams = UpdateStruct(ModemParams, Modem);
  end
end

% derived quantities
code  = readldpccode(ModemParams.Codec);
H     = ldpccode2matrix(code);
if ModemParams.InfoLen~=(size(H,2)-size(H,1)); error('wrong length'); end;
ModemParams.K              = size(H,2)-size(H,1);
ModemParams.Hu             = H(:,1:ModemParams.K);

ModemParams.Fs            = ModemParams.SymbolRate*ModemParams.P;

ModemParams.pilotindex    = ModemParams.pilot.pos-1+(1:length(ModemParams.pilot.bits));

ModemParams.CodedSymLenR  = (ModemParams.CodewordLen+length(ModemParams.pilot.bits))/log2(length(ModemParams.ModulationMap));
ModemParams.TxRrcCoefs    = rcosfir(ModemParams.RrcAlpha, ModemParams.RrcNsym, ModemParams.P, 1, 'sqrt');
ModemParams.RxRrcCoefs    = ModemParams.TxRrcCoefs;

mexQAMModulator(0, 'init', length(ModemParams.pilot.bits)/log2(length(ModemParams.ModulationMap)), ModemParams.ModulationMap);
ModemParams.pilot.syms    = mexQAMModulator(0,'modulate', ModemParams.pilot.bits);

% ModemParams.TguardSam     = ModemParams.SlotLenSam-((ModemParams.CodedSymLenR+ModemParams.RrcNsym)*ModemParams.P);
ModemParams.Tguard        = ModemParams.TguardSam/ModemParams.Fs;
ModemParams.MaxTd         = ModemParams.Tguard;
ModemParams.SlotLenSam    = ModemParams.TguardSam + 2*ModemParams.RrcNsym*ModemParams.P+1 + (ModemParams.CodedSymLenR-1)*ModemParams.P;

%% Acquisition Parameters
STAcqParams.SSR           = 7;
STAcqParams.fsin          = 1/(2*(1/ModemParams.SymbolRate)) * (1+ModemParams.RrcAlpha/2);
STAcqParams.ToneRampSym   = 2; % duration of ramp at the beginning and end of the tone in symbols
STAcqParams.FilterTone    = 0;
STAcqParams.CancelTone    = 1;


STAcqParams.taumin        = 0;
STAcqParams.taumax        = ModemParams.Tguard;
STAcqParams.f_min         = -2;
STAcqParams.f_max         = +2;
STAcqParams.fr_min        = -60;
STAcqParams.fr_max        = -10;

STAcqParams.pshapeName    = 'rrc';
STAcqParams.rolloff       = ModemParams.RrcAlpha;
T                         = (1/ModemParams.SymbolRate);
Ts                        = T/ModemParams.P;
STAcqParams.pshape        = @(t) rrcpulse(t,ModemParams.RrcAlpha,T,ModemParams.P) .* (abs(t) <= ModemParams.RrcNsym*ModemParams.P*Ts);

STAcqParams.timingEstimatorName = 'TdEst';
STAcqParams.phaseEstimatorName  = 'PhEst';

STAcqParams.Lsin          = ModemParams.CodedSymLenR;
STAcqParams.amplsin       = [];

STAcqParams.M             = length(ModemParams.ModulationMap);
STAcqParams.T             = 1/ModemParams.SymbolRate;
STAcqParams.OS            = ModemParams.P;
STAcqParams.Tpshape       = 2*ModemParams.RrcNsym/ModemParams.SymbolRate;
STAcqParams.P             = ModemParams.pilot.pos-1 + (1:(length(ModemParams.pilot.bits)/log2(length(ModemParams.ModulationMap))));
tmp                       = ones(1,ModemParams.CodedSymLenR);
tmp(STAcqParams.P)        = 0;
STAcqParams.D             = find(tmp);
STAcqParams.p             = ModemParams.pilot.syms;

STAcqParams.pTime     = STAcqParams.p;
STAcqParams.pPhase    = STAcqParams.p;
STAcqParams.DTime     = STAcqParams.D;
STAcqParams.PTime     = STAcqParams.P;
STAcqParams.DPhase    = STAcqParams.D;
STAcqParams.PPhase    = STAcqParams.P;
STAcqParams.searchOS  = 2;
STAcqParams.Nfft      = 8192;

% update default structure with user supplied values
if nargin>=2
  if ~isempty(STAcq)
    STAcqParams = UpdateStruct(STAcqParams, STAcq);
  end
end

% derived quantities
STAcqParams.SSR           = STAcqParams.SSR+10*log10(ModemParams.P);
STAcqParams.amplsin       = sqrt(1/(10^(STAcqParams.SSR/10)));
