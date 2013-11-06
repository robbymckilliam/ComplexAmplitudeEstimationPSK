% scramble data
%---------------------------------------------------------------------
% 
% ScrambledData = scramble_bits(ScramblerState, DataIn)
%
% Inputs
%   ScramblerState : Apply this as initial state of scrambler
%                    Data bits are scrambled using the 
%                    polynomial x^7 + x^4 + 1 
%   DataIn         : Data to be scrambled
%
% Outputs
%   ScrambledData  : Scrambled bits
% Comments:
%   This function is used to scramble and un-scramble data bits
%   note that scramble(State,scramble(State,Data)) = Data
% 
%---------------------------------------------------------------------
% Author:  David Haley
% project: ASRP
%---------------------------------------------------------------------
% Copyright 2012  
% Institute for Telecommunications Research
% University of South Australia
%---------------------------------------------------------------------
% 
%---------------------------------------------------------------------

function ScrambledData = scramble_bits(ScramblerState, DataIn)

ShiftReg = ScramblerState;
ScrambledData = zeros(length(DataIn),1);

for BitIndex = 1:length(DataIn)
  
  ScrambleBit = xor(ShiftReg(4),ShiftReg(7));
  ScrambledData(BitIndex) = xor(DataIn(BitIndex),ScrambleBit);
  ShiftReg(2:end) = ShiftReg(1:end-1);
  ShiftReg(1) = ScrambleBit;
end




