function c  = conv_fft(a, b)
%compute a convolution using the fft.  The output is identical to:
%conv(a,b,`valid')  when the length(a) >= length(b),  and
%conv(b,a,`valid') when length(a) =< length(b).
%
% (c) 2012 ITR-UniSA
% Authors: Robby McKilliam
% Created: 3 April 2012
    
    P = numel(a);
    Q = numel(b);
    L = P + Q - 1;
    M = min([P,Q]);

    c = ifft(fft(a,L) .* fft(b,L));
    c = c(M:end-M+1);

end