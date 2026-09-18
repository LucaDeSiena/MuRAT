function [energyRatioBodyCoda_i,energyRatioCodaNoise_i]   =...
   Murat_body(bodyWindow,startNoise,srate_i,sp_i,...
   cursorPick_i,cursorCodaStart_i,cursorCodaEnd_i)
% function [energyRatioBodyCoda_i,energyRatioCodaNoise_i] =...
% Murat_body(bodyWindow,startNoise,srate_i,cursorPick_i,...
% cursorCodaStart_i,cursorCodaEnd_i)
%
% CREATES the body-to-coda energy ratios and uncertainties necessary for
%   the CN method.
%
% Input parameters:
%    bodyWindow:                window for body wave measurements
%    startNoise:                start of noise window
%    srate_i:                   sampling rate
%    sp_i:                      filtered amplitude
%    cursorPick_i:              picking along the trace
%    cursorCodaStart_i:         start of coda window
%    cursorCodaEnd_i:           end of coda window    
%
% Output parameters:
%    energyRatioBodyCoda_i:     energy ratio between body and coda waves
%    energyRatioCodaNoise_i:    energy ratio between coda waves and noise

% Redefine windows with sampling for noise and direct waves
intSamples              =   floor(bodyWindow*srate_i);
cursorNoise0            =	floor(startNoise*srate_i);
cursor0Noise1           =   cursorNoise0 + intSamples;
cursorPS1               =   floor(cursorPick_i + intSamples-1);
lengthSamples           =   length(cursorCodaStart_i:cursorCodaEnd_i);

% Extract windows (vectorized)
spPS                    =   sp_i(cursorPick_i:cursorPS1,:);
spNoise                 =   sp_i(cursorNoise0:cursor0Noise1,:);
spCoda                  =   sp_i(cursorCodaStart_i:cursorCodaEnd_i,:);

% Integrate using trapz with correct spacing
dt                      =   1 / srate_i;
spampPS                 =   trapz((0:dt:dt*(size(spPS,1)-1)), spPS); 
spampNoise              =   trapz((0:dt:dt*(size(spNoise,1)-1)), spNoise);


% Normalize coda by relative window length
differentWindows        =   lengthSamples/intSamples;
spampCoda               =   ...
    trapz((0:dt:dt*(size(spCoda,1)-1)), spCoda)./ differentWindows;

energyRatioBodyCoda_i   =   spampPS ./  spampCoda;
energyRatioCodaNoise_i  =   spampCoda ./ spampNoise;
    
end