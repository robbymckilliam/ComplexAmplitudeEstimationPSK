function [rrc] = rrcpulse(t,rolloff,T,OS)

% root raised cosine pulse
%
% INPUT:
%   t       - vector of time instants at which pulse is evaluated
%   rolloff - rolloff factor
%   T       - symbol period
%   OS      - oversampling factor

Ts = T/OS;    % sample period
fc = 1/(2*T); % cut off frequency
fs = 1/Ts;    % sampling frequency

% implementation based on MATLAB's firrcos

% use logical indexing:
tol = 2*eps;
idx1 = (abs(t) < tol);
% idx2 = (abs(abs(8*rolloff*fc*t) - 1.0) < sqrt(eps));
idx2 = (abs(abs(8*rolloff*fc*t) - 1.0) < tol);
idx3 = ~idx1 & ~idx2;

rrc = zeros(size(t));
rrc(idx1) = -sqrt(2.*fc) ./ (pi.*fs) .* (pi.*(rolloff-1) - 4.*rolloff );
rrc(idx2) = sqrt(2.*fc) ./ (2.*pi.*fs) ...
		* (    pi.*(rolloff+1)  .* sin(pi.*(rolloff+1)./(4.*rolloff)) ...
		- 4.*rolloff .* sin(pi.*(rolloff-1)./(4.*rolloff)) ...
		+ pi.*(rolloff-1) .* cos(pi.*(rolloff-1)./(4.*rolloff)) ...
		);
rrc(idx3) = -4.*rolloff./fs .* ( cos((1+rolloff).*2.*pi.*fc.*t(idx3)) + ...
        sin((1-rolloff).*2.*pi.*fc.*t(idx3)) ./ (8.*rolloff.*fc.*t(idx3)) ) ...
        ./ (pi .* sqrt(1./(2.*fc)) .* ((8.*rolloff.*fc.*t(idx3)).^2 - 1));
    
% scale pulse
rrc = rrc * sqrt(2.*fc);

end