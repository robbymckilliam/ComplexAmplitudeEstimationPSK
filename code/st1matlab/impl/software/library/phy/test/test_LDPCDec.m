function success = test_LDPCDec()

success = false;

addpath ../bin

rng('default');
rng(0);

% display the module under test
fullpath = which('mexLDPCDec');
if isempty(fullpath)
  disp('could not find mexLDPCDec');
  return;
end
fprintf('testing %s\n', fullpath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LDGM ENCODER TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read LDGM code definition
try
  mexLDPCDec(0, 'readdecoder', 'LDGM168.dec');
catch
  disp('could not read decoder definition');
  return;
end

% check if successfully loaded
if ~mexLDPCDec(0, 'decloaded')
  disp('decoder should be loaded by now');
  return;
end

% get dimensions
[N,M] = mexLDPCDec(0, 'dimensions');
if (N~=588) || (M~=420)
  disp('wrong dimensions');
  return;  
end

% generate random data
info    = randi(2,1,N-M)-1;

% read code definition in Matlab and compute parity
code    = readldpccode('LDGM168.dec');
H       = ldpccode2matrix(code);
Hu      = H(:,1:168);
parity2 = rem(info*Hu',2);

% test double data type
% encode
parity = mexLDPCDec(0, 'encodeLDGM', double(info));
if length(parity)~=M
  disp('LDGM encoder returned wrong length using doubles');
  return;
end
if ~strcmp(class(parity),'double')
  disp('LDGM encoder returned wrong type using doubles');
  return;  
end
if ~isequal(sort(unique(parity)), [0 1])
  disp('LDGM parity contains values other than 0 and 1 using doubles');
  return;
end

% compute and verify parity
if ~isequal(parity, parity2)
  disp('LDGM wrong encoding using doubles');
  return;
end

% test single data type
% encode
parity = mexLDPCDec(0, 'encodeLDGM', single(info));
if length(parity)~=M
  disp('LDGM encoder returned wrong length using singles');
  return;
end
if ~strcmp(class(parity),'single')
  disp('LDGM encoder returned wrong type using singles');
  return;  
end
if ~isequal(sort(unique(parity)), [0 1])
  disp('LDGM parity contains values other than 0 and 1 using singles');
  return;
end

% compute and verify parity
if ~isequal(parity, parity2)
  disp('LDGM wrong encoding using singles');
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RA ENCODER TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read RA code definition
try
  mexLDPCDec(0, 'readdecoder', 'RA168.dec');
catch
  disp('could not read decoder definition');
  return;
end

% check if successfully loaded
if ~mexLDPCDec(0, 'decloaded')
  disp('decoder should be loaded by now');
  return;
end

% get dimensions
[N,M] = mexLDPCDec(0, 'dimensions');
if (N~=336) || (M~=168)
  disp('wrong dimensions');
  return;  
end

% generate random data
info    = randi(2,1,N-M)-1;

% read code definition in Matlab and compute parity
code    = readldpccode('RA168.dec');
H       = ldpccode2matrix(code);
Hu      = H(:,1:168);
parity2 = rem(cumsum(rem(info*Hu',2)),2);

% test double data type
% encode
parity = mexLDPCDec(0, 'encodeRA', double(info));
if length(parity)~=M
  disp('RA encoder returned wrong length using doubles');
  return;
end
if ~strcmp(class(parity),'double')
  disp('RA encoder returned wrong type using doubles');
  return;  
end
if ~isequal(sort(unique(parity)), [0 1])
  disp('RA parity contains values other than 0 and 1 using doubles');
  return;
end

% compute and verify parity
if ~isequal(parity, parity2)
  disp('RA wrong encoding using doubles');
  return;
end

% test single data type
% encode
parity = mexLDPCDec(0, 'encodeRA', single(info));
if length(parity)~=M
  disp('RA encoder returned wrong length using singles');
  return;
end
if ~strcmp(class(parity),'single')
  disp('RA encoder returned wrong type using singles');
  return;  
end
if ~isequal(sort(unique(parity)), [0 1])
  disp('RA parity contains values other than 0 and 1 using singles');
  return;
end

% compute and verify parity
if ~isequal(parity, parity2)
  disp('RA wrong encoding using singles');
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECODER TEST (USES RA ENCODER)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read RA code definition (already tested before)
mexLDPCDec(0, 'readdecoder', 'RA168.dec');

% get dimensions
[N,M]   = mexLDPCDec(0, 'dimensions');

% generate random data
info    = randi(2,1,N-M)-1;

% encode
parity  = mexLDPCDec(0, 'encodeRA', double(info));

% modulate to BPSK and add Gaussian noise
sigma   = 0.7;
maxit   = 100;
x       = [info parity];
tx      = 1-2*x;
y       = tx + sigma*randn(1,N);
Lch     = 2*y/sigma^2;

% decode with doubles
[Lapp, it] = mexLDPCDec(0, 'decode', double(Lch), maxit);
if length(Lapp)~=N
  disp('decoder returned wrong length using doubles');
  return;
end
if ~strcmp(class(Lapp),'double')
  disp('decoder returned wrong type using doubles');
  return;  
end
if sum(isnan(Lapp))>0
  disp('detected NaN after decoding using doubles');
  return;
end
if sum(isinf(Lapp))>0
  disp('detected Inf after decoding using doubles');
  return;
end
if sum((Lapp<0)~=x)>0
  disp('decoder was unable to decode using doubles');
  return;
end
if it==maxit
  disp('decoder decoded successful but did not stop using doubles');
  return;
end

% decode with singles
[Lapp, it] = mexLDPCDec(0, 'decode', single(Lch), maxit);
if length(Lapp)~=N
  disp('decoder returned wrong length using singles');
  return;
end
if ~strcmp(class(Lapp),'single')
  disp('decoder returned wrong type using singles');
  return;  
end
if sum(isnan(Lapp))>0
  disp('detected NaN after decoding using singles');
  return;
end
if sum(isinf(Lapp))>0
  disp('detected Inf after decoding using singles');
  return;
end
if sum((Lapp<0)~=x)>0
  disp('decoder was unable to decode using singles');
  return;
end
if it==maxit
  disp('decoder decoded successful but did not stop using singles');
  return;
end

disp('all tests passed');
success = true;

end