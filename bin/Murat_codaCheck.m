function [tCoda_i, cursorCodaStart_i, cursorCodaEnd_i]   =...
    Murat_codaCheck(originTime_i,pktime_i,srate_i,tCm,tWm,tempis,...
    peakDelay_i,peakDelayMethod)

% function [tCoda_i, cursorCodaStart_i, cursorCodaEnd_i]   =...
%     Murat_codaCheck(originTime_i,pktime_i,srate_i,tCm,tWm,tempis,...
%     peakDelay_i,peakDelayMethod)
%
% SETS the correct lapse time. In case it is not defined you are
%   measuring envelopes from the peak of the direct wave assuming that
%   the entire waveform is diffusive.
%
% Input parameters:
%    originTime_i:          origin time in seconds
%    pktime_i:              piked time in seconds
%    srate_i:               sampling rate
%    tCm:                   coda starting time
%    tWm:                   length of coda window
%    tempis:                time vector from seismogram    
%    peakDelay_i:           peak delay for coda lapse time    
%    peakDelayMethod:       choise between constant, peak, or ravel time    
%
% Output parameters:
%    tCoda_i:               coda starting time after check in seconds
%    cursorCodaStart_i:     coda starting time after check on trace
%    cursorCodaEnd_i:       coda end time after check on trace

%Define envelope duration
t00                 =   tempis(1);
nSamples            =   numel(tempis);

% compute reference offsets once
pk_rel              =   pktime_i - originTime_i;

% Method peak is generally valid for active seismics
switch peakDelayMethod
    case 'Peak'
        tCoda_i     =   pk_rel + peakDelay_i;
    case 'Constant'
        tCoda_i     =   tCm;
    case 'Travel'
        tCoda_i     =   originTime_i + tCm * pk_rel;
    otherwise
        error("Unknown peak-delay method")
end

% Define the indexes along the seismogram

cursorCodaStart_i   =   floor((tCoda_i + originTime_i - t00) * srate_i -1);

endSample0          =   cursorCodaStart_i + floor(tWm * srate_i - 1);
cursorCodaEnd_i     =   min(nSamples, endSample0);

end
