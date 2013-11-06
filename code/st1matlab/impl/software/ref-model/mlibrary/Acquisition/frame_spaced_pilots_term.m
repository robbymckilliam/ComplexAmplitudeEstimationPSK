% Assemble frame of data with one or more pilot segments
%--------------------------------------------------------------------------
% frame_spaced_pilots_term.m: This function assembles a frame with one 
% or more pilot segment(s) at user-specified positions and pilot bits.
%
% [frame,idx_pilots,idx_data] = frame_spaced_pilots_term(pilot,...
%                                     data,modparam,AISmode,softmode)
% 
%
% Inputs:
%   pilot(n)    - array of structs with the following fields for each pilot
%                 segment:
%       .pos        - position index where pilot segment starts
%       .bits       - vector of pilot bits {0,1}
%   data        - vector of binary data {0,1}.
%                 Optionally, data can the string 'dummy' to indicate that
%                 the data portion of the frame is to be filled with dummy
%                 data (all ones). This can be useful when the data portion
%                 is not of interest, e.g. when constructing a local
%                 receiver copy of the pilot signal for acquisition
%                 purposes. When data='dummy', no data bits are appended
%                 after the last pilot segment.
%   modparam    - struct holding the modulation parameters
%       .T          - symbol period
%       .BT         - GMSK bandwidth/time product
%       .P          - oversampling rate
%       .L          - GMSK memory length
%   AISmode     - [OPTIONAL] true to run in AIS mode, false otherwise.
%                 Default value is false.
%   softmode    - [OPTIONAL] true if using LLRs instead of bits, false
%                 otherwise. Default value is false.
% Outputs:
%   frame       - assembled frame with pilot, data and termination segments
%   idx_pilots  - vector of pilot position indices
%   idx_data    - vector of data position indices
% Comments:
% Assemble a frame with one or more pilot segment(s) at user-specified
% positions and pilot bits. The remaining bit positions are then filled up
% with data and the trellis is terminated before each pilot segment. For
% details on how the termination is implemented, see gmsk_termtrellis.m.
%
% By default, this function operates in non-AIS mode. When run in AIS-mode,
% the input argument pilot must contain the AIS training and start flag as
% first pilot segment. Additional pilot segments in the data portion of the
% AIS frame are possible. Before returning the frame, the AIS preamble is
% removed as it will later be added by AISframe.m. Furthermore, a reserved
% zero-bit is assumed before the trellis termination of any additional
% pilot segment. This reserved bit ensures that the termination bits will
% not be affected by AIS bit stuffing.
%
% By default, this function expects binary {0,1} data and pilots as input.
% However, when softmode is set to true, the input parameters data and
% pilot(n).bits can be used to pass log-likelihood ratios (LLR) instead of
% {0,1}. In this case, the LLRs of the trellis termination bits are set to
% zero, i.e. they assumed to be unknown apriori.
%--------------------------------------------------------------------------
% Author:  A. Pollock, ITR
% project: ASRP
%--------------------------------------------------------------------------
% Copyright 2011 : 
% Institute for Telecommunications Research
% University of South Australia
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
function [frame,idx_pilots,idx_data] = frame_spaced_pilots_term(...
                                        pilot,data,modparam,AISmode,softmode)

if (nargin < 4)
    AISmode = false;
end
if (nargin < 5)
    softmode = false;
end

% number of pilot segments
N_pilots = length(pilot);
% number of termination bits
if (mod(modparam.L,2) == 0)
    N_term_bits = modparam.L + 1;
else
    N_term_bits = modparam.L;
end
if AISmode
    % an additional reserved bit to avoid bitstuffing is required
    N_term_bits = N_term_bits + 1;
end
if strcmp(modparam.ModScheme,'qpsk'),N_term_bits=0;end;
% compute pilot start and end positions
pilot_start_pos = zeros(1,N_pilots);
N_pilot_bits = zeros(1,N_pilots);
pilot_bits = [];
for n = 1:N_pilots
    pilot_start_pos(n) = pilot(n).pos;
    N_pilot_bits(n) = length(pilot(n).bits);
    pilot_bits = [pilot_bits pilot(n).bits];
end
pilot_end_pos = pilot_start_pos + N_pilot_bits - 1;

% compute data start and end positions
data_start_pos = pilot_end_pos(1:end-1) + 1;
% check if data can be inserted before first pilot symbol
if (pilot_start_pos(1) > 1)
    if (pilot_start_pos(1) > N_term_bits+1)
        data_start_pos = [1 data_start_pos];
    else
        error('Insufficient space to insert data and termination bits before the first pilot bit.')
    end
end
data_end_pos = pilot_start_pos - N_term_bits - 1;
% exclude non-positive indices
data_end_pos = data_end_pos(data_end_pos > 0);
if ~isempty(data_end_pos)
    N_data_bits = data_end_pos - data_start_pos + 1;
else
    N_data_bits = [];
end
if ~strcmp(data,'dummy')
    % check if any data bits remain; append to frame if any
    if (sum(N_data_bits) < length(data))
        data_start_pos(end+1) = pilot_end_pos(end)+1;
        N_data_bits(end+1) = length(data) - sum(N_data_bits);
        data_end_pos(end+1) = data_start_pos(end) + N_data_bits(end) - 1;
    end
end
% number of data segments
N_data = length(data_start_pos);
% number of trellis terminations
N_term = sum(data_end_pos < pilot_end_pos(end));

%% assemble frame
frame = nan(1,sum(N_pilot_bits)+sum(N_data_bits)+N_term*N_term_bits);

if strcmp(data,'dummy')
    % insert dummy data
    if ~AISmode
        % assume all-ones data
        data = ones(1,sum(N_data_bits));
    else
        % assume all-ones data, which still needs to be anti-bit-stuffed
        data_tmp = ones(1,sum(N_data_bits));
    end
else
    if (length(data) < sum(N_data_bits))
        error('Insufficient data to fill bit positions between pilot segments.');
    end
end

if AISmode
    % define AIS training and flag
    training = repmat([0 1], 1, 12);
    flag     = [0 1 1 1 1 1 1 0];
    % check that first pilot is AIS preamble
    if ~isequal(pilot(1).bits,[training flag]) || (pilot(1).pos ~= 1)
        error('In AIS mode pilot(1) must be the AIS training and start flag.');
    end
    % check that data segments between pilots and termination bits are
    % integer multiples of bytes
    if any(rem(N_data_bits,8)) || any(rem(N_pilot_bits(2:end)+N_term_bits,8))
        error('In AIS mode all data segments and termination+pilot segments must be integer multiples of bytes.');
    end
    
    if ~strcmp(data,'dummy')
        % apply AIS byte flipping to data bits
        data = ByteFlipping(data);
    else
        % flip ones in dummy data to avoid bitstuffing in AISframe.m
        data = AntiBitStuffing(data_tmp);
    end
end

% insert pilot bits
idx_pilots = [];
for n = 1:N_pilots
    idx_pilots = [idx_pilots, pilot_start_pos(n):pilot_end_pos(n)];
end
frame(idx_pilots) = pilot_bits;

% insert data bits
idx_data = [];
for n = 1:N_data
    idx_data = [idx_data, data_start_pos(n):data_end_pos(n)];
end
 frame(idx_data) = data;

% trellis termination 
if ~softmode && ~strcmp(modparam.ModScheme,'qpsk')
    % compute and insert trellis termination bits
    % term_start_pos = data_end_pos(data_end_pos < pilot_end_pos(end)) + 1;
    for pos = data_end_pos(data_end_pos < pilot_end_pos(end))
        if ~AISmode
            % standard binary mapping {0,1}->{+1,-1}
            frame_mapped = 1-2*frame(1:pos);
            [term_bits_mapped] = gmsk_termtrellis(frame_mapped,modparam.L);
            % inverse mapping to obtain {0,1} data
            frame_with_term_bits = (1-[frame_mapped, term_bits_mapped])/2;
            frame(pos+[1:N_term_bits]) = frame_with_term_bits(pos+[1:N_term_bits]);
        else
            % append reserved 0 bit to avoid bitstuffing of termination bits
            frame_tmp = [frame(1:pos), 0];
            % AIS NRZI encoding
            frame_nrzi = AISNRZI(frame_tmp);
            [term_bits_nrzi] = gmsk_termtrellis(frame_nrzi,modparam.L);
            % inverse mapping to obtain {0,1} data
            frame_with_term_bits = AISNRZIinv([frame_nrzi, term_bits_nrzi]);
            frame(pos+[1:N_term_bits]) = frame_with_term_bits(pos+[1:N_term_bits]);
        end
    end
else
    % assume no prior knowledge about termination bits
    frame(isnan(frame)) = 0;
end

if AISmode
    % remove preamble since this gets added by AISframe.m
    frame      = frame(length([training flag])+1:end);
    idx_pilots = idx_pilots(length([training flag])+1:end) - length([training flag]);
    idx_data   = idx_data - length([training flag]);
    % apply inverse AIS byte flipping to frame
    frame = ByteFlipping(frame);
    
%     warning('idx_pilots refers to the AIS byte flipped version of frame.')
end

end