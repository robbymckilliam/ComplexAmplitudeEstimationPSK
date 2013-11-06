function [UserParams, TxPackets] = st1transmitter(UserParams, ModemParams, STAcqParams)

%-------------------------------------------------------------%
%               Encode/Interleave Data Source
%-------------------------------------------------------------%

UserParams.Data = encoder(scramble_bits(ones(1,8), UserParams.Info)', UserParams);             
UserParams.Data = interleave(UserParams.Data, UserParams.S);
%-------------------------------------------------------------%
%               Insert Pilot Information Bits
%-------------------------------------------------------------%
if UserParams.N_pilots                
  UserParams.Data = frame_spaced_pilots_term(UserParams.pilot, UserParams.Data, ModemParams); 
end            
%-------------------------------------------------------------%
%           Transmit Modulation 
%-------------------------------------------------------------%
Syms      = mexQAMModulator(0, 'modulate', UserParams.Data);
TxPackets = pulseshapefilter(Syms, ModemParams);

if ModemParams.AddTone %~ModemParams.PerfectChan
  TxPackets = addtone(TxPackets, UserParams, ModemParams, STAcqParams);
end

UserParams.TxSyms = Syms; %only used for perfect channel SNR est

end
