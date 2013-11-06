function [phihat] = phaseest_da(RxSyms,acqparams)

% extract received symbols corresponding to pilots
RxSyms = RxSyms(acqparams.P);

C = RxSyms .* conj(acqparams.p);
C = sum(C)/length(acqparams.p);

phihat = angle(C);