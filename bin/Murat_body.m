function [energyRatioBodyCoda_i,energyRatioCodaNoise_i]...
                        =   Murat_body(bodyWindow,startNoise,srate_i,...
                    sp_i,cursorPick_i,cursorCodaStart_i,cursorCodaEnd_i)
% function [energyRatioBodyCoda_i,energyRatioCodaNoise_i] =...
% Murat_body(Murat,srate_i,sp_i,cursorPick_i,cursorCodaStart_i,cursorCodaEnd_i)
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

int                     =   bodyWindow*srate_i;
lc                      =   length(cursorCodaStart_i:cursorCodaEnd_i);
cursor0                 =	floor(startNoise*srate_i);
cursor0_1               =   floor(cursor0 + int);
cursor2                 =   floor(cursorPick_i + int-1);

spamp                   =...
    trapz(sp_i(cursorPick_i:cursor2,:))/int;

spampn                  =...
    trapz(sp_i(cursor0:cursor0_1,:))/int;

spampc                  =...
    trapz(sp_i(cursorCodaStart_i:cursorCodaEnd_i,:))/lc;

energyRatioBodyCoda_i   =   spamp./spampc;
energyRatioCodaNoise_i  =   spampc./spampn;
    
end