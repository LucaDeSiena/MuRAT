function cursorPeakDelay_i      =...
    Murat_peakDelayCheck(tempis,cursorPick_i,maxtpde,srate_i)
% function cursorPeakDelay_i      =...
%     Murat_peakDelayCheck(tempis,cursorPick_i,maxtpde,srate_i)
%
% CHECKS that the peak delay is inside the trace and warns if not
%
% Input parameters:
%    tempis:            time vector of seismogram
%    cursorPick_i:      picking on the trace
%    maxtpde:           maximum time to look for peak
%    srate_i:           sampling rate
%    
% Output parameters:
%    cursorPeakDelay_i: location of peak delay window

cursorPeakDelayMax              =   floor(cursorPick_i+maxtpde*srate_i);
lengthTempis                    =   length(tempis);

%Measures for peak delay - compute shorter window if waveform is cut
if cursorPeakDelayMax > lengthTempis
    
    cursorPeakDelay_i           =   lengthTempis; 
    msg                         =...
        'Length of seismogram shorter than peak-delay window!';
    warning(msg);
    
else
    
    cursorPeakDelay_i           =   cursorPeakDelayMax;
    
end          


end
