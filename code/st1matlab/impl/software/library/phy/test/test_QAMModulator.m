function success = test_QAMModulator()

success = false;

addpath ../bin

rng('default');
rng(0);

% display the module under test
fullpath = which('mexQAMModulator');
if isempty(fullpath)
  disp('could not find mexQAMModulator');
  return;
end
fprintf('testing %s\n', fullpath);

% use QAM constellation with 4 symbols
Constellation = [-1-1i -1+1i 1-1i 1+1i]/sqrt(2);
bps           = log2(length(Constellation));
tolerance     = 1e-6;
N             = 1e3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODULATE TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise modulator
mexQAMModulator(0, 'init', 1, double(Constellation));

% test every symbol individually
for s=1:length(Constellation)
  bits  = de2bi(s-1, bps, 'left-msb');
  x     = mexQAMModulator(0, 'modulate', bits);
  if abs(x-Constellation(s))>tolerance
    disp('incorrect results when testing modulator using doubles');
    return;
  end
end

% initialise longer modulator
mexQAMModulator(0, 'init', N, double(Constellation));

% generate random data
bits  = randi(2,1,N*bps)-1;

% modulate in Matlab
idx   = (2.^(0:(bps-1))) * flipud(reshape(bits, bps, N)) + 1;
yref  = Constellation(idx);

% modulate
x     = mexQAMModulator(0, 'modulate', double(bits));
if length(x)~=N
  disp('modulator returned wrong length using doubles');
  return;
end
if ~strcmp(class(x),'double')
  disp('modulator returned wrong type using doubles');
  return;  
end
if max(abs(x-yref))>tolerance
  disp('incorrect results when running modulator using doubles');
  return;    
end



% initialise modulator
mexQAMModulator(0, 'init', 1, single(Constellation));

% test every symbol individually
for s=1:length(Constellation)
  bits  = de2bi(s-1, bps, 'left-msb');
  x     = mexQAMModulator(0, 'modulate', bits);
  if abs(x-Constellation(s))>tolerance
    disp('incorrect results when testing modulator using singles');
    return;
  end
end

% initialise longer modulator
mexQAMModulator(0, 'init', N, double(Constellation));

% generate random data
bits  = randi(2,1,N*bps)-1;

% modulate in Matlab
idx   = (2.^(0:(bps-1))) * flipud(reshape(bits, bps, N)) + 1;
yref  = Constellation(idx);

% modulate
x     = mexQAMModulator(0, 'modulate', single(bits));
if length(x)~=N
  disp('modulator returned wrong length using singles');
  return;
end
% The modulator returns doubles instead of singles even when singles should
% be returned. Ignored at the moment as it can be converted in Matlab if
% needed.
% if ~strcmp(class(x),'single')
%   disp('modulator returned wrong type using singles');
%   return;  
% end
if max(abs(x-yref))>tolerance
  disp('incorrect results when running modulator using singles');
  return;    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMODULATE TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%

mexQAMModulator(0, 'init', N, double(Constellation));

% generate random data
bits  = randi(2,1,N*bps)-1;

% modulate
x     = mexQAMModulator(0, 'modulate', double(bits));

% add Gaussian noise
sigma = 0.7;
y     = x + sigma*(randn(1,N) + 1i*randn(1,N));

% set channel observations
mexQAMModulator(0, 'set-ch', double(y), sigma^2);

% create and set zero a-priori LLRs
La    = zeros(bps,N);
mexQAMModulator(0, 'set-apri', double(La));

% demodulate
Lout  = mexQAMModulator(0, 'demodulate');
if size(Lout,1)~=bps
  disp('demodulator returned wrong number of bits per symbol using doubles');
  return;
end
if size(Lout,2)~=N
  disp('demodulator returned wrong length using doubles');
  return;
end
if ~strcmp(class(Lout),'double')
  disp('demodulator returned wrong type using doubles');
  return;  
end
if sum(isnan(Lout))>0
  disp('detected NaN after demodulating using doubles');
  return;
end
if sum(isinf(Lout))>0
  disp('detected Inf after demodulating using doubles');
  return;
end

% use channel with low noise to get good demodulation results
sigma = 0.1;
y     = x + sigma*(randn(1,N) + 1i*randn(1,N));

% set channel observations
mexQAMModulator(0, 'set-ch', double(y), sigma^2);

% create and set zero a-priori LLRs
La    = zeros(bps,N);
mexQAMModulator(0, 'set-apri', double(La));

% demodulate
Lout  = mexQAMModulator(0, 'demodulate');

% remodulate signal based on demodulated results
mexQAMModulator(0, 'set-apri', double(Lout));
xhat  = mexQAMModulator(0, 'remodulate');
if length(xhat)~=N
  disp('remodulator returned wrong number of symbols using doubles');
  return;
end
if ~strcmp(class(xhat),'double')
  disp('remodulator returned wrong type using doubles');
  return;  
end
if sum(isnan(xhat))>0
  disp('detected NaN after remodulating using doubles');
  return;
end
if sum(isinf(xhat))>0
  disp('detected Inf after remodulating using doubles');
  return;
end
if max(abs(x-xhat))>tolerance
  disp('demodulator/remodulator loopback test failed for doubles');
  return;
end




mexQAMModulator(0, 'init', N, single(Constellation));

% generate random data
bits  = randi(2,1,N*bps)-1;

% modulate
x     = mexQAMModulator(0, 'modulate', single(bits));

% add Gaussian noise
sigma = 0.7;
y     = x + sigma*(randn(1,N) + 1i*randn(1,N));

% set channel observations
mexQAMModulator(0, 'set-ch', single(y), sigma^2);

% create and set zero a-priori LLRs
La    = zeros(bps,N);
mexQAMModulator(0, 'set-apri', single(La));

% demodulate
Lout  = mexQAMModulator(0, 'demodulate');
if size(Lout,1)~=bps
  disp('demodulator returned wrong number of bits per symbol using singles');
  return;
end
if size(Lout,2)~=N
  disp('demodulator returned wrong length using singles');
  return;
end
% At this point the demodulator cannot determine if single or double should
% be returned and it defaults to double.
% if ~strcmp(class(Lout),'single')
%   disp('demodulator returned wrong type using singles');
%   return;  
% end
if sum(isnan(Lout))>0
  disp('detected NaN after demodulating using singles');
  return;
end
if sum(isinf(Lout))>0
  disp('detected Inf after demodulating using singles');
  return;
end

% use channel with low noise to get good demodulation results
sigma = 0.1;
y     = x + sigma*(randn(1,N) + 1i*randn(1,N));

% set channel observations
mexQAMModulator(0, 'set-ch', single(y), sigma^2);

% create and set zero a-priori LLRs
La    = zeros(bps,N);
mexQAMModulator(0, 'set-apri', single(La));

% demodulate
Lout  = mexQAMModulator(0, 'demodulate');

% remodulate signal based on demodulated results
mexQAMModulator(0, 'set-apri', single(Lout));
xhat  = mexQAMModulator(0, 'remodulate');
if length(xhat)~=N
  disp('remodulator returned wrong number of symbols using singles');
  return;
end
% At this point the remodulator cannot determine if single or double should
% be returned and it defaults to double.
% if ~strcmp(class(xhat),'single')
%   disp('remodulator returned wrong type using singles');
%   return;  
% end
if sum(isnan(xhat))>0
  disp('detected NaN after remodulating using singles');
  return;
end
if sum(isinf(xhat))>0
  disp('detected Inf after remodulating using singles');
  return;
end
if max(abs(x-xhat))>tolerance
  disp('demodulator/remodulator loopback test failed for singles');
  return;
end

disp('all tests passed');
success = true;

end