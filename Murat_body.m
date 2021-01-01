function [energyRatioBodyCoda_i,energyRatioCodaNoise_i]...
                        = Murat_body(Murat,srate_i,hsp_i,...
                        cursorPick_i,cursorCodaStart_i,cursorCodaEnd_i)
%CREATES the body-wave-to-coda energy ratios and uncertainties

% Noise and direct-wave amplitude
lc                      =   length(cursorCodaStart_i:cursorCodaEnd_i);
cursor0                 =   floor(Murat.input.startNoise*srate_i);
int                     =   Murat.input.bodyWindow*srate_i;% samples, P/S
intn                    =   int;% samples, noise
cursor0_1               =   floor(cursor0+intn);
cursor2                 =   floor(cursorPick_i + int-1); %end-sample

spamp                   =...
    trapz(abs(hsp_i(cursorPick_i:cursor2)))/int;% the direct

spampn                  =...
    trapz(abs(hsp_i(cursor0:cursor0_1)))/intn;% the noise

spampc                  =...
    trapz(abs(hsp_i(cursorCodaStart_i:cursorCodaEnd_i)))/lc;% the coda

energyRatioBodyCoda_i   =   spamp/spampc;%spectral ratios signal-coda
energyRatioCodaNoise_i  =   spampc/spampn;%spectral ratios coda-noise
