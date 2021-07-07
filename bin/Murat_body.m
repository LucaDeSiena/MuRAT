function [energyRatioBodyCoda_i,energyRatioCodaNoise_i]...
                        =   Murat_body(Murat,srate_i,sp_i,...
         cursorPick_i,cursorCodaStart_i,cursorCodaEnd_i)

% CREATES the body-to-coda energy ratios and uncertainties
int                     =   Murat.input.bodyWindow*srate_i;
lc                      =   length(cursorCodaStart_i:cursorCodaEnd_i);
cursor0                 =	floor(Murat.input.startNoise*srate_i);
cursor0_1               =   floor(cursor0+int);
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