function [phihat] = phaseest_viterbi(RxSyms,M)

F = 1; % F could also be selected as a function of RxSyms
C = F .* (RxSyms./abs(RxSyms)).^M;
C = sum(C)/length(RxSyms);

phihat = angle(C)/M;