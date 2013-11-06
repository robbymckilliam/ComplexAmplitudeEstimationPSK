% Updates elements of one struct to another
%--------------------------------------------------------------------------
% out = UpdateStruct(in, update)
%
% Inputs:
%   in       struct that will be updated
%   update   struct which holds the values to be updated
%
% Outputs:  
%   out      updated struct
% 
%--------------------------------------------------------------------------
% Author: Gottfried Lechner
% Project: ASRP
%--------------------------------------------------------------------------
% Copyright 2013
% Institute for Telecommunications Research
% University of South Australia
%--------------------------------------------------------------------------

function out = UpdateStruct(in, update)

out   = in;
names = fieldnames(update);

for i=1:length(names)
  out.(names{i}) = update.(names{i});
end

end
