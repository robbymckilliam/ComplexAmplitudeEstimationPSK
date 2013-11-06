function H = ldpccode2matrix(code)
%ldpccode2matrix
%
%   converts LDPC code struct to sparse parity-check matrix
%
%   H = ldpccode2matrix(code)
%
%   H       ... sparse parity-check matrix
%
%   code    ... LDPC code struct
%
%   (c) 2010 Gottfried Lechner / ITR-UniSA.
%   gottfried.lechner@unisa.edu.au

N = length(code.vardegree);
M = length(code.chkdegree);
I = length(code.interleaver);

c = zeros(N, max(code.vardegree));

m = zeros(1,I);
o = 1;
for i=1:N
    m(o:o+code.vardegree(i)-1) = i;
    o = o+code.vardegree(i);
end

m = m(code.interleaver);

o = 1;
for i=1:M
    var = m(o:o+code.chkdegree(i)-1);
    for j=1:size(var,2)
        c(var(j),find(c(var(j),:)==0, 1)) = i;
    end
    o = o+code.chkdegree(i);
end

rows = max(max(c));
cols = size(c,1);

H = spalloc(rows,cols,length(nonzeros(c)));

for i=1:cols
   H(nonzeros(c(i,:)),i) = 1;
end

end
