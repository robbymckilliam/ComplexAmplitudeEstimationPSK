% This function simply removes an index of pilots from the data
%--------------------------------------------------------------------------
% RemovePilots.m - This function simply removes an index of pilots 
% from the data. It also check if pilots are expected in data portion 
% of packet and if not, is transparent and returns output identical to
% input.
% 
% [ dataout] = RemovePilots( datain,userparam,modparam )
% Inputs:
%   datain      - Data sequence including pilots
%   userparam   - User Information containing pilot index
%   modparam    - Modem Information (e.g. AIS mode etc...)
% Outputs:
%   dataout     - Data sequence with pilots removed
% Comments:
%      
%--------------------------------------------------------------------------
% Author:   M. Lavenant, ITR
% Project : ASRP
%--------------------------------------------------------------------------
% Copyright 2011 : 
% Institute for Telecommunications Research
% University of South Australia
%--------------------------------------------------------------------------
% Contains ITR Background IP
%--------------------------------------------------------------------------
function [ dataout] = RemovePilots( datain,~,modparam )

dataout=datain;
if modparam.AISFormat
    if modparam.N_pilots >1 %more than just default start train/flag
        startoffset=modparam.TrainL+modparam.FlagL;
        dataout( modparam.pilotindex(startoffset+1:end)-startoffset)=[];
    end
else
    if modparam.N_pilots >0
        dataout( modparam.pilotindex)=[];
    end
end

