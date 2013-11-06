function SoftDec = ST1ULacquisitionWrapper(SoftDec, rx, acqparams, ModemParams)

% % set acquisition parameters
% acqparams.pshapeName          = 'rrc';
% acqparams.rolloff             = ModemParams.RrcAlpha;
% acqparams.timingEstimatorName = 'TdEst';
% acqparams.phaseEstimatorName  = 'PhEst';
% acqparams.fr_max              = ModemParams.MaxFd;
% acqparams.fr_min              = 0;
% 
% % initialise acquisition (should be moved outside)
% ST1ULacquisition_init(acqparams);

% call actual acquisition module
[mfoutput, fhat, frhat, amplhat, tauhat, nhat] = ST1ULacquisition(rx, acqparams);

% change referencepoint from the middle of the packet to the first sample
% in the slot
deltaT = -tauhat-length(mfoutput)/2/ModemParams.SymbolRate - (length(ModemParams.TxRrcCoefs)-1)/2/ModemParams.Fs;

% transform estimates based on new reference point
[phihat,fhat,frhat] = polynomialPhaseTransform(angle(amplhat), fhat, frhat, deltaT);
phihat = mod(phihat, 2*pi);

%-----------------------------------------------------------------%
%       Transfer channel offset estimates to SD structure
%-----------------------------------------------------------------%
if tauhat<0 % check validity of estimates (to avoid indexing errors)
  tauhat = 0;
elseif tauhat>acqparams.taumax
  tauhat = acqparams.taumax;
end  
SoftDec.Fo     = fhat/ModemParams.Fs;
SoftDec.Fd     = frhat/ModemParams.Fs;
SoftDec.Td     = tauhat*ModemParams.Fs;
SoftDec.PhOff  = phihat;

% scaling of signal and noise power
scaling             = ModemParams.P*ModemParams.Fs.^2;
SoftDec.SigPowEst   = abs(amplhat).^2 * scaling; 
SoftDec.NoiseVar    = nhat*2 * scaling;

SoftDec.Acquired    = 1;
SoftDec.RxSyms      = mfoutput;
end
