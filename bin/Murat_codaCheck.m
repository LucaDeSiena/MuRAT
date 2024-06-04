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
t00                                             =   tempis(1);
lengthTempis                                    =   length(tempis);

% Method peak is generally valid for active seismics
if isequal(peakDelayMethod,'Peak')
    
    tCoda_i                                     =...
        (pktime_i-originTime_i)+peakDelay_i;

% Method constant is the most used in small local tomography 
elseif isequal(peakDelayMethod,'Constant')
    
    tCoda_i                                     =   tCm;

% Method travel is the standard at the regional scale
elseif isequal(peakDelayMethod,'Travel')
    
    tCoda_i                                     =...
        originTime_i+ tCm*(pktime_i-originTime_i);

end

% Define the indexes along the seismogram
cursorCodaStart_i                               =...
    floor((originTime_i - t00 + tCoda_i) * srate_i - 1);

cursorCodaEnd_i                                 =...
    floor(cursorCodaStart_i + tWm * srate_i - 1);

% Some seismograms are not long enough so you consider a shorter window
if cursorCodaEnd_i > lengthTempis
    cursorCodaEnd_i                             =   lengthTempis;
end

end
