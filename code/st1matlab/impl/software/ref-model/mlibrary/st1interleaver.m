function S = st1interleaver()

% parity bits through block interleaver
P = 1:(8*29);
P = reshape(P, 8, 29);
P = reshape(P', 1, 8*29);

% interleave with information bits
S = reshape([1:232 ; P+232], 1, 464);

end
