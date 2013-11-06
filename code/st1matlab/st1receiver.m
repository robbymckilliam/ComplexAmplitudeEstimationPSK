function [SoftDec, NoiseHypothesis] = st1receiver(RxSamples, ModemParams, STAcqParams, SoftDec)

% various debug plot settings - active only one at a time
debugplot   = 0;  % shows spectrum and time/frequency position of packets
scatterplot = 0;  % scatterdiagram
timeplot    = 0;  % time-domain plot comparing received and reconstructed waveforms

%% check parameters
% we should not filter and cancel the tone
if STAcqParams.FilterTone && STAcqParams.CancelTone
  warning('the tone gets filterd and cancelled, i.e., removed twice');
end

%% pre-compute parameters for tone cancellation
if STAcqParams.CancelTone
  % Compute tone attenuation caused by the root-raised cosine
  % receive filter.
  RRCAttenuation  = sqrt(1/2*(1+cos(pi/ModemParams.SymbolRate/ModemParams.RrcAlpha*(STAcqParams.fsin - (1-ModemParams.RrcAlpha)*ModemParams.SymbolRate/2))));

  % Signal to sinusoid ratio (symbol level)
  SSRSym          = 10^(STAcqParams.SSR/20)/sqrt(ModemParams.P);

  % relative amplitude of the tone after matched filtering
  shiftampl       = RRCAttenuation/SSRSym;

  % set up symbol times
  tsym            = (0:(ModemParams.CodedSymLenR-1))/ModemParams.SymbolRate;
  
  % precompute the effect of the tone on the symbols
  ToneInSymbols   = shiftampl*exp(+2i*pi*STAcqParams.fsin*tsym);
  
  STAcqParams.ToneInSymbols = ToneInSymbols;
end


%% Initialise Soft Decoders
if nargin>3
  if isempty(SoftDec)
    clear SoftDec;
  end
end

if exist('SoftDec', 'var')
  startsd = length(SoftDec)+1;
else
  startsd = 1;
end

for sd=startsd:ModemParams.SoftDecoders
  SoftDec(sd) = initialiseSoftDec(ModemParams.CodedSymLenR, length(RxSamples)); %#ok<*AGROW>
end
  

%% Multi user detector
NoiseHypothesis = RxSamples;
AllDone         = 0;

for m=1:ModemParams.MUDIter  
  %% reset the received signals of all soft decoders
  [SoftDec.RxSyms] = deal([]);

  %% Acquisition
  % skip acquisition if exisiting SoftDecs have been provided and we are in
  % the first iteration
  if ~((m==1) && (nargin >= 4))
    % if we have an empty soft decoder then we try to acquire a new signal
    idx = find([SoftDec.Acquired]==0,1);
    if idx
      % only acquire if all running soft decoders cannot make any further progress
      StillRunning = 0;
      for sd=1:ModemParams.SoftDecoders
        if(SoftDec(sd).Acquired) && (~SoftDec(sd).success)
          StillRunning = StillRunning || ((SoftDec(sd).DecStat-SoftDec(sd).DecStatP)>0.10); % 0.05 -> 0.10
        end
      end

      if ~StillRunning
        % acquire a new user and put it in a free soft decoder
        SoftDec(idx)     = ST1ULacquisitionWrapper(SoftDec(idx), NoiseHypothesis, STAcqParams, ModemParams);      
        SoftDec(idx).TTL = ModemParams.TTL;

        % TODO: check if we should really accept this user; ideally using a
        % confidence interval of the timing estimate and/or tone
      end
    end
  end
  
  %% Soft Decoder
  ICSamples = NoiseHypothesis;
  
  % sort/randomise active soft decoders
  idx     = find([SoftDec.Acquired] & ~[SoftDec.finish]);
  [~,ix]  = sort([SoftDec(idx).SigPowEst]./[SoftDec(idx).NoiseVar], 'descend');
%   ix      = randperm(length(idx));
  idx     = idx(ix);
  
  for sic=1:length(idx)
    sd = idx(sic);

    % save decoder state
    SoftDec(sd).DecStatP  = SoftDec(sd).DecStat;

    for refineiter=1:ModemParams.RefineIter
      if isempty(SoftDec(sd).RxSyms)
        % if the user was not acquired in this iteration then we need to rebuild the signal
        if ~isempty(SoftDec(sd).SampHatp)
          % we cancelled this user already -> add the contribution back
          ICSamples = ICSamples + SoftDec(sd).SampHatp;
        end

        % compensate for channel effects and normalise
        ModemParams.ChanType  = 'chancomp';
        DataComp              = satchan(ICSamples, SoftDec(sd), ModemParams, SoftDec(sd), [],  STAcqParams);
        SoftDec(sd).RxSyms    = DataComp/sqrt(SoftDec(sd).SigPowEst);
      end

      % debug plot of scatter diagram
      if scatterplot; hold off; scatter(real(SoftDec(sd).RxSyms), imag(SoftDec(sd).RxSyms), 'rx'); axis(1.5*[-1 1 -1 1]); axis square; grid on; title(sprintf('Decoder %i', sd)); end %#ok<*UNRCH>

      % cancellation of tone in symbol domain
      % !!! It seems that inside the initial acquisition we also cancel the
      % tone if this flag is set. Hence, we cancel the tone twice. This
      % shouldn't be a problem in subsequent refinement iterations as the
      % tone will be added when rebuilding the signal but is certainly
      % sub-optimal in the first iteration. -> FIX IT 2013/05/02 Andre,
      % Gottfried
      if STAcqParams.CancelTone
        % energy before tone cancellation
        Ebefore         = mean(abs(SoftDec(sd).RxSyms.^2));

        % just subtract the precomputed influence of the tone
        SoftDec(sd).RxSyms  = SoftDec(sd).RxSyms - ToneInSymbols;

        % energy after tone cancellation
        Eafter          = mean(abs(SoftDec(sd).RxSyms.^2));

        % renormalise the energy in the signal
        % this seems to be required because we end up with a negative SNR
        % otherwise - needs to be investigated further
        SoftDec(sd).RxSyms  = SoftDec(sd).RxSyms * sqrt(Ebefore/Eafter);

        if scatterplot; hold on; scatter(real(SoftDec(sd).RxSyms), imag(SoftDec(sd).RxSyms), 'bx'); end
      end

      % debug plot of scatter diagram
      if scatterplot; pause; end        

      %-------------------------------------------------------------%
      % Estimate noise&interference variance
      %-------------------------------------------------------------%
      % Assume the matched filtered received signal for the user of
      % interest is normalised by its transmit power, i.e., Px = 1.
      % Measure the received power - everything in excess of Px is
      % noise&interference.
      % TODO: fix
%       SNRest  = 1/(mean(abs(SoftDec(sd).RxSyms).^2)-1);
      SNRest = SoftDec(sd).SigPowEst/SoftDec(sd).NoiseVar;
      if SNRest<0
%           fprintf('negative SNR estimated (%.2f) for decoder %i\n', SNRest, sd);
        SNRest = 2;
      end
      SNRest  = min(SNRest, 100);
      sigma2  = 1/SNRest;

      % run the demodulator/decoder for QAM
      if ~SoftDec(sd).success
        SoftDec(sd) = SingleUserRXQAM(SoftDec(sd), sigma2, ModemParams, sd);
      end

      % decide if we should drop this decoder
      if (refineiter==ModemParams.RefineIter) && (~SoftDec(sd).success)
        SoftDec(sd).TTL = SoftDec(sd).TTL - 1;
        if SoftDec(sd).TTL==0
%           fprintf('removed decoder %i (expired TTL)\n', sd);
          SoftDec(sd) = initialiseSoftDec(ModemParams.CodedSymLenR, length(RxSamples));
          break;
        end
      end

      %-------------------------------------------------------------%
      % Refine estimates
      %-------------------------------------------------------------%
      
      if ~strcmpi(STAcqParams.phaseEstimatorName,'disable')
          STAcqParamsRefine = STAcqParams;
          STAcqParamsRefine.P = zeros(1,length(STAcqParams.P)+length(STAcqParams.D));
          STAcqParamsRefine.P(STAcqParams.P) = STAcqParams.P;
          STAcqParamsRefine.P(STAcqParams.D) = STAcqParams.D;
          STAcqParamsRefine.p = zeros(1,length(STAcqParams.P)+length(STAcqParams.D));
          STAcqParamsRefine.p = SoftDec(sd).SoftSyms; % soft symbols for pilots + data
          PhOff = phaseest_da(SoftDec(sd).RxSyms,STAcqParamsRefine);
          maxcorr = abs(PhOff);
          SoftDec(sd).PhOff = PhOff;
      else
          maxcorr = 0;
      end
      
% %       maxcorr = 0;
%       [SoftDec(sd), maxcorr] = ST1ULacquisitionRefinePeriodoWrapper(ICSamples, SoftDec(sd), ModemParams, STAcqParams);
% %       SoftDec(sd) = ST1ULacquisitionRefineNoiseVar(ICSamples, SoftDec(sd), ModemParams, STAcqParams);


      if (SoftDec(sd).success) && ((maxcorr<0.01) || (ModemParams.SoftDecoders==1)) % no further refinement possible if the only decoder was successful...
        SoftDec(sd).finish = 1;
        break;
      end
      SoftDec(sd).RxSyms = [];

      %-------------------------------------------------------------%
      % Rebuild signal
      %-------------------------------------------------------------%
      if SoftDec(sd).Acquired
        TxSampHat = pulseshapefilter(SoftDec(sd).SoftSyms, ModemParams);
        TxSampHat = addtone(TxSampHat, [], ModemParams, STAcqParams);

        ModemParams.ChanType  = 'chanadd';
        SoftDec(sd).SampHat   = satchan(TxSampHat, SoftDec(sd), ModemParams, SoftDec(sd));      

        ICSamples = ICSamples - SoftDec(sd).SampHat;
        SoftDec(sd).SampHatp  = SoftDec(sd).SampHat;

        % compare the received and the reconstructed waveforms
        if timeplot; duration = 1000; subplot(1,2,1); hold off; plot(real(RxSamples(1:duration)), 'b'); grid on; hold on; plot(real(SoftDec(sd).SampHat(1:duration)),'r--'); subplot(1,2,2); hold off; plot(real(RxSamples((end-duration):end)), 'b'); grid on; hold on; plot(real(SoftDec(sd).SampHat((end-duration):end)),'r--'); pause; end
      end
    end
  end
  
  %% update noise hypothesis
  NoiseHypothesis = RxSamples;
  for sd=1:ModemParams.SoftDecoders
    if SoftDec(sd).Acquired
      SoftDec(sd).SampHatp  = SoftDec(sd).SampHat;
      NoiseHypothesis       = NoiseHypothesis - SoftDec(sd).SampHat;
      
%       % DEBUG: filter tone
%       ts = (0:(length(NoiseHypothesis)-1))/ModemParams.Fs;
%       NoiseHypothesis = NoiseHypothesis .* exp(-2i*pi*(STAcqParams.fsin+SoftDec(sd).Fo*ModemParams.Fs)*ts);
%       NoiseHypothesis = dcblockST1UL(NoiseHypothesis);
%       NoiseHypothesis = NoiseHypothesis .* exp(+2i*pi*(STAcqParams.fsin+SoftDec(sd).Fo*ModemParams.Fs)*ts);
%       % DEBUG END
    end
  end
  
  % TODO: check the power spectral density of noisehypothesis and decide
  % whether to continue
  
  if debugplot; st1receiver_debugplot; end
  
%   save(sprintf('debug%03i.mat', m));
    
  %% cleanup duplicates
  for sd=1:(ModemParams.SoftDecoders-1)
    if SoftDec(sd).success
      for j=(sd+1):ModemParams.SoftDecoders
        if (SoftDec(j).success) && (sum(SoftDec(sd).RxMsgDec~=SoftDec(j).RxMsgDec)==0)
%           fprintf('removed decoder %i (same message as decoder %i)\n', j, sd);
          SoftDec(j) = initialiseSoftDec(ModemParams.CodedSymLenR, length(RxSamples));
        end
      end
    end
  end

  if AllDone
    break;
  end

  % all decoders have decoded - nothing more to do
  % but do one more iteration to clean up the residual signal
  if all([SoftDec.finish])
    AllDone = 1;
  end  
end

% update noise hypothesis
NoiseHypothesis = RxSamples;
for sd=1:ModemParams.SoftDecoders
  if SoftDec(sd).finish
    NoiseHypothesis = NoiseHypothesis - SoftDec(sd).SampHat;
  end
end  

end
