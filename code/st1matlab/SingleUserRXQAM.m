function SoftDec = SingleUserRXQAM2(SoftDec, sigma2, ModemParams, SoftDecID)

% set channel observations for demodulator
mexQAMModulator(0, 'set-ch', SoftDec.RxSyms, sigma2);
mexQAMModulator(0, 'set-apri', zeros(log2(length(ModemParams.ModulationMap)), length(SoftDec.RxSyms)));

%---------------------------------------------------------------------%
%               DEMODULATE SIGNAL
%---------------------------------------------------------------------%
Lext  = mexQAMModulator(0, 'demodulate');
Ladec = Lext(:)';

%---------------------------------------------------------------------%
%               PILOTS REMOVAL (transparent if no pilots)
%---------------------------------------------------------------------%
Ladec = RemovePilots(Ladec, [], ModemParams);
%---------------------------------------------------------------------%
%               De-Interleave (transparent if no interleaver)
%---------------------------------------------------------------------%
Ladec = deinterleave(Ladec, ModemParams.S);
%---------------------------------------------------------------------%
%              DECODE SIGNAL
%---------------------------------------------------------------------%
[Lapp,it]       = mexLDPCDec(SoftDecID, 'decode', Ladec, ModemParams.maxit);
Lu              = Lapp(1:ModemParams.InfoLen);
SoftDec.success = (it<ModemParams.maxit);
if SoftDec.success
  Ledec = 100*sign(Lapp);
  Lapp  = 100*sign(Lapp);
else
  Ledec = Lapp - Ladec;
end
SoftDec.DecStat = mean(tanh(abs(Lapp)/2));
if sum(isnan(Ledec))
  disp('STOP - NaNs detected in Ledec');
  keyboard
end

% use either Lapp or Ledec for rebuild
Lout = Ledec;
% Lout = Lapp;

%---------------------------------------------------------------------%
%               INTERLEAVER (transparent if no interleaver)
%---------------------------------------------------------------------%
Lout = interleave(Lout, ModemParams.S);
%---------------------------------------------------------------------%
%               PILOT INSERTION (if required)
%---------------------------------------------------------------------%
if ModemParams.N_pilots>0
  dummy       = ModemParams.pilot; %replace data pilots with LLRs
%   meanRel     = mean(abs(Lout));
%   meanRel     = 1; %%% WRONG: SHOULD BE SET TO MAXLLR (e.g. 100)
  meanRel = 100; %%% Should be correct. Will test if it has a big impact. (Gottfried, Andre, 2013/09/03)
  
  for n=1:ModemParams.N_pilots
    dummy(n).bits = meanRel*(1-2*ModemParams.pilot(n).bits);
  end
  Lout = frame_spaced_pilots_term(dummy,Lout,ModemParams,0,1);
end
%---------------------------------------------------------------------%
%               remodulate signal using ext and app
%---------------------------------------------------------------------%
mexQAMModulator(0, 'set-apri', (reshape(Lout, log2(length(ModemParams.ModulationMap)), [])));
SoftDec.SoftSyms = mexQAMModulator(0, 'remodulate');  

%---------------------------------------------------------------------%
%               Make Hard Decision on LLRs
%---------------------------------------------------------------------%
if SoftDec.success
  SoftDec.RxMsgDec    = scramble_bits(ones(1,8), Lu<0)';
else
  SoftDec.RxMsgDec    = zeros(1,length(Lu));
end

end