function st1uplink_su_monte_carlo(ResFile,EbNodB_vec_in,framesPerSnr_in,minNumFrameErr_in,ModemParams_in,STAcqParams_in)

% Script to run ST1 uplink in single-user mode

%% simulation parameters

EbNodB_vec          = [0:0.5:10];
framesPerSnr        = 1e6;
minNumFrameErr      = 1e2;
numFramesBeforeSave = 5e3;

%% non-default ST1 system parameters

ModemParams.P             = 4;

ModemParams.AddTone       = 0;
ModemParams.pilot.pos     = 1;
% ModemParams.pilot.bits    = [0 0 1 1 0 0 1 1 0 0 0 0 1 1 0 1 1 0 0 1 1 0 0 1 0 1 1 0]; % ST1 pilot (14 symbols)
% ModemParams.pilot.bits    = [0 0 1 1 0 0 1 1 0 0 0 0 1 1 0 1 1 0 0 1 1 0 0 1 0 1 1 0, 1 1]; % 15-symbol pilot
% ModemParams.pilot.bits    = [0 0 1 1 0 0 1 1 0 0 0 0 1 1 0 1 1 0 0 1]; % 10-symbol pilot
ModemParams.pilot.bits    = [0 0 1 1 0 0 1 1 0 0]; % 5-symbol pilot

ModemParams.SoftDecoders  = 1;
ModemParams.MUDIter       = 1;
% without turbo sync:
%{%
ModemParams.maxit         = 100; % max LDPC decoder iterations
ModemParams.RefineIter    = 1;   % turbo sync iterations
%}
% with turbo sync:
%{
ModemParams.maxit         = 25; % max LDPC decoder iterations
ModemParams.RefineIter    = 4;  % turbo sync iterations
%}

%% non-default ST1 acquisition parameters

STAcqParams.timingEstimatorName   = 'disable';
% STAcqParams.timingEstimatorName = 'TdEst';

% STAcqParams.phaseEstimatorName    = 'disable';
% STAcqParams.phaseEstimatorName    = 'data-aided';
STAcqParams.phaseEstimatorName    = 'Mackenthun';

STAcqParams.SSR           = Inf; % i.e. no tone
STAcqParams.CancelTone    = 0;

STAcqParams.taumin        = 0;
STAcqParams.taumax        = 0;
STAcqParams.f_min         = 0;
STAcqParams.f_max         = 0;
STAcqParams.fr_min        = 0;
STAcqParams.fr_max        = 0;
STAcqParams.fr_min_refine = 0;
STAcqParams.fr_max_refine = 0;

%% set up results file name (with timestamp)

if (nargin >= 1)
  ResFile = [ResFile '_' datestr(now,'yyyymmddTHHMMSS')];
  if exist('./results','dir')
    ResFile = ['results/' ResFile];
  else
    error('Directory ''./results'' does not exist.')
  end
else
  ResFile = [];
end

if ~isdeployed
  setpath
end

%% overwrite default parameters with user-supplied parameters

if nargin>=2
  EbNodB_vec = EbNodB_vec_in;
end
if nargin>=3
  framesPerSnr = framesPerSnr_in;
end
if nargin>=4
  minNumFrameErr = minNumFrameErr_in;
end

% update structures with user-supplied values
if nargin>=5
  if ~isempty(ModemParams_in)
    ModemParams = UpdateStruct(ModemParams, ModemParams_in);
  end
end
if nargin>=6
  if ~isempty(STAcqParams_in)
    STAcqParams = UpdateStruct(STAcqParams, STAcqParams_in);
  end
end
clear ModemParams_in STAcqParams_in
% save non-default parameters
if ~isempty(ResFile)
  save(ResFile);
end
[ModemParams, STAcqParams] = overwriteDefaultParams(ModemParams, STAcqParams);

%% initialise acquisition, modulator and decoder(s)
ST1ULacquisition_init(STAcqParams);
mexQAMModulator(0, 'init', ModemParams.CodedSymLenR, ModemParams.ModulationMap);
for k=1:ModemParams.SoftDecoders
  mexLDPCDec(k, 'readdecoder', ModemParams.Codec);
end

%% channel parameter offsets
% time delay is uniformly distributed between Td_min and Td_max [sec]
Td_min = 0;
Td_max = 0;

% set to true to activate uniformly distributed phase offsets between -/+pi
% rand_phase_flag = false;
rand_phase_flag = true;

% frequency offset is uniformly distributed between -/+Fo_max [Hz]
Fo_max = 0;

% Doppler rate is uniformly distributed between Fd_min and Fd_max [Hz/s]
Fd_min = 0;
Fd_max = 0;

if ~isempty(ResFile)
  save(ResFile, ...
    'Td_min','Td_max','rand_phase_flag','Fo_max','Fd_min','Fd_max', ...
    '-append');
end

%% Monte-Carlo simulation loop

% initialise global random number stream based on current time
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

% Eb/No to Es/No conversion
% code rate and codeword length relative to total number of bits
CodeRate        = ModemParams.InfoLen/ModemParams.CodewordLen;
CodewordLenRel  = ModemParams.CodewordLen / (ModemParams.CodewordLen + length(ModemParams.pilot.bits));

EbNodB_incUW = EbNodB_vec - 10*log10(CodewordLenRel);
EsNodB_vec   = EbNodB_vec + 10*log10(CodeRate*log2(length(ModemParams.ModulationMap)));

numFramesSimulated  = zeros(1,length(EbNodB_vec));
numBitErr           = zeros(1,length(EbNodB_vec));
numFrameErr         = zeros(1,length(EbNodB_vec));

mseerrtau   = zeros(1,length(EbNodB_vec));
mseerrfreq  = zeros(1,length(EbNodB_vec));
mseerrfrate = zeros(1,length(EbNodB_vec));
mseerrampl  = zeros(1,length(EbNodB_vec));
mseerrphase = zeros(1,length(EbNodB_vec));

printstr = ['Eb/No [dB]', ...
  '\t', ' FER ', ...
  '\t', ' RMSEs:', ...
  '\t', ' Td [%%]', ...
  '\t\t', ' Fo [Hz]',...
  '\t', ' Fr [Hz/s]',...
  '\t',' Ph [rad]', ...
  '\t',' Ampl'];
numCharFprint = fprintf([printstr, '\n']);
fprintf([repmat('=',1,length(printstr)+(length(printstr)-numCharFprint)*3), '\n']);

for EbNoind = 1:length(EbNodB_vec)
  
  tic
  for frameind = 1:framesPerSnr
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   TRANSMITTER                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % generate random information bits
    ModemParams.Info = randi(2,[1,ModemParams.InfoLen])-1;
    
    % create transmit packet
    [~, Tx]          = st1transmitter(ModemParams, ModemParams, STAcqParams);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                      CHANNEL                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % generate random channel offsets
    PhOff   = pi*(2*rand-1) * rand_phase_flag;
    Td      = Td_min + rand * (Td_max-Td_min);  % time offset
    Fo      = (2*rand-1) * Fo_max;              % carrier offset
    Fd      = Fd_min + rand * (Fd_max-Fd_min);  % Doppler rate
    
    % fixed unit channel gain
    Ampl    = 1;
    
    ChanParams.PhOff        = PhOff;
    ChanParams.Td           = Td*ModemParams.Fs;
    ChanParams.Fo           = Fo/ModemParams.Fs;
    ChanParams.Fd           = Fd/ModemParams.Fs;
    ChanParams.FdDir        = 1;
    ChanParams.SigPowEst    = Ampl*10^(EsNodB_vec(EbNoind)/10); % Note that ST1 framework assumes unit noise variance
    
    if strcmpi(STAcqParams.phaseEstimatorName,'disable') || strcmpi(STAcqParams.phaseEstimatorName,'data-aided')
      % if complex gain estimator disabled, genie aid signal power estimate
      STAcqParams.SigPowEst = ChanParams.SigPowEst;
    end
    
    % simulate channel
    ModemParamsChan = ModemParams;
    ModemParamsChan.ChanType = 'chanadd';
    Rx = satchan(Tx, ChanParams, ModemParamsChan, [], []);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     RECEIVER                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [SoftDec] = st1receiver(Rx, ModemParams, STAcqParams);
    
    % replace above line by the following code to simulate uncoded BER
    %{
        % assume zero time offset, but adjust search window:
        %   The timing estimator assumes the first sample at t=Ts and the first
        %   symbol at t=T=OS*Ts, i.e. it does not take the leading Tx filter tail
        %   into account. Therefore, we need to shift the search window by the
        %   leading Tx filter tail minus (OS-1) sample periods.
        STAcqParams.Tpshape = 2*ModemParams.RrcNsym/ModemParams.SymbolRate;
        STAcqParams.T   = 1/ModemParams.SymbolRate;
        STAcqParams.Ts  = 1/ModemParams.Fs;
        STAcqParams.OS = ModemParams.P;
        STAcqParams.pshape  = @(t) rrcpulse(t,0.5,STAcqParams.T,ModemParams.P) .* (abs(t) <= ModemParams.RrcNsym*ModemParams.P*STAcqParams.Ts);
        tauhatMF = 0 + STAcqParams.Tpshape/2 - (ModemParams.P-1)*STAcqParams.T/ModemParams.P;
        mfoutput = mf(Rx, tauhatMF, STAcqParams);
        mexQAMModulator(0, 'init', length(mfoutput), ModemParams.ModulationMap);
        mexQAMModulator(0, 'set-ch', mfoutput, 1);
        Lext  = mexQAMModulator(0, 'demodulate');
        Ladec = Lext(:)';
        ModemParams.S = st1interleaver();
        Ladec(ModemParams.pilotindex) = [];
        Ladec = deinterleave(Ladec, ModemParams.S);
        % extract info bits
        Ladec = Ladec(1:ModemParams.InfoLen);
        SoftDec.RxMsgDec = scramble_bits(ones(1,8), Ladec<0)';
        
        SoftDec.PhOff = 0;
        SoftDec.Td = 0;
        SoftDec.Fo = 0;
        SoftDec.Fd = 0;
        SoftDec.SigPowEst = ChanParams.SigPowEst;
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     EVALUATION                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % count bit and frame errors
    numBitErrCur = sum(SoftDec.RxMsgDec ~= ModemParams.Info);
    if (numBitErrCur > 0)
      numBitErr(EbNoind)   = numBitErr(EbNoind) + numBitErrCur;
      numFrameErr(EbNoind) = numFrameErr(EbNoind) + 1;
      if (numFrameErr(EbNoind) >= minNumFrameErr)
        numFramesSimulated(EbNoind) = frameind;
        break
      end
    end
    
    % estimation errors
    PhOffHat    = SoftDec.PhOff;
    TdHat       = SoftDec.Td/ModemParams.Fs;
    FoHat       = SoftDec.Fo*ModemParams.Fs;
    FdHat       = SoftDec.Fd*ModemParams.Fs;
    AmplHat     = SoftDec.SigPowEst/ChanParams.SigPowEst;
    
    mseerrtau(EbNoind)   = mseerrtau(EbNoind)   + (Td - TdHat)^2;
    mseerrfreq(EbNoind)  = mseerrfreq(EbNoind)  + (Fo - FoHat)^2;
    mseerrfrate(EbNoind) = mseerrfrate(EbNoind) + (Fd - FdHat)^2;
    mseerrampl(EbNoind)  = mseerrampl(EbNoind)  + (Ampl - AmplHat)^2;
    mseerrphase(EbNoind) = mseerrphase(EbNoind) + (angle(exp(1i*PhOff)*exp(-1i*PhOffHat)))^2;
    
    % save intermediate results
    if ~mod(frameind,numFramesBeforeSave) && ~isempty(ResFile)
      save(ResFile, ...
        'EbNoind','frameind', ...
        'mseerrtau','mseerrfreq','mseerrfrate','mseerrampl','mseerrphase', ...
        'numFramesSimulated','numBitErr','numFrameErr', ...
        '-append');
    end
    
  end
  numFramesSimulated(EbNoind) = frameind;
  
  elapsedTime = toc;
  fprintf([num2str(EbNodB_vec(EbNoind)),...
    '\t\t ', num2str(numFrameErr(EbNoind)/numFramesSimulated(EbNoind),2),...
    '\t\t ', num2str(100*sqrt(mseerrtau(EbNoind)/numFramesSimulated(EbNoind))*ModemParams.SymbolRate,2),...
    '\t\t ', num2str(sqrt(mseerrfreq(EbNoind)/numFramesSimulated(EbNoind)),2), ...
    '\t\t ', num2str(sqrt(mseerrfrate(EbNoind)/numFramesSimulated(EbNoind)),2), ...
    '\t\t ', num2str(sqrt(mseerrphase(EbNoind)/numFramesSimulated(EbNoind)),2), ...
    '\t\t ', num2str(sqrt(mseerrampl(EbNoind)/numFramesSimulated(EbNoind)),2), ...
    '\t (', num2str(elapsedTime,'%.1f'), ' sec)', '\n']);
  
  % save intermediate results
  if ~isempty(ResFile)
    save(ResFile, ...
      'EbNoind', ...
      'mseerrtau','mseerrfreq','mseerrfrate','mseerrampl','mseerrphase', ...
      'numFramesSimulated','numBitErr','numFrameErr', ...
      '-append');
  end
  
end

BER = numBitErr./numFramesSimulated/ModemParams.InfoLen;
FER = numFrameErr./numFramesSimulated;
if any(numFrameErr < minNumFrameErr)
  warning('The minimum number of frame errors has not been reached at all SNRs.')
end
% BERuncodTheory = berawgn(EbNodB_vec+10*log10(CodeRate),'psk',length(ModemParams.ModulationMap),'nondiff');

mseerrtau       = mseerrtau./numFramesSimulated;
relrmserrtau    = 100*sqrt(mseerrtau)*ModemParams.SymbolRate; % relative RMS estimation error in percent
mseerrfreq      = mseerrfreq./numFramesSimulated;
mseerrfrate     = mseerrfrate./numFramesSimulated;
mseerrampl      = mseerrampl./numFramesSimulated;
mseerrphase     = mseerrphase./numFramesSimulated;

% save final results
if ~isempty(ResFile)
  save(ResFile, ...
    'mseerrtau','relrmserrtau','mseerrfreq','mseerrfrate','mseerrampl','mseerrphase', ...
    'numFramesSimulated','numBitErr','numFrameErr','BER','FER', ...
    '-append');
end

end