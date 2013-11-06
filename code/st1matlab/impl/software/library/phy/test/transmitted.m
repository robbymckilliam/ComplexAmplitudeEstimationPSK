function out = transmitted(t, s, params)
%generate an amplited modulated transmitted signal
    g = params.g;
    T = params.T;
    P = params.P;
    D = params.D;
    csum = zeros(size(t));
    for k = 1:length(D)
        csum = csum + s(D(k))*g(t - D(k)*T);
    end
    for k = 1:length(P)
        csum = csum + s(P(k))*g(t - P(k)*T);
    end
    out = csum;
       
end