function TxSymbols = pulseshapefilter(Syms, ModemParams)

%-------------------------------------------------------------%
%   Pulse Shape Filter (RRC)
%-------------------------------------------------------------%
TxSymbols = upfirdn(Syms,ModemParams.TxRrcCoefs,ModemParams.P,1);


% now the tail of the filter is included to match the hardware
% implementation (5.2.2013 - Andre, Gottfried)
%--------------------------------------------------------%
%----------- Truncate Rrc Filter end Tail----------------%
% Ltail     = length(TxSymbols)-(length(Syms)*ModemParams.P);
% Lh        = length(ModemParams.TxRrcCoefs);
% Lend      = Ltail-(Lh-1)/2;
% TxSymbols(end-Lend+1:end)=[];%chop end tail only

end
