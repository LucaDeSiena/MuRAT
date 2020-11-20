function [theoreticalTime_i, originTime_i, tC_i, cursorPick_i,...
    cursorPeakDelay_i, cursorCodaStart_i, cursorCodaEnd_i] =...
    Murat_times(Murat,tempis,pktime_i,sst_i,srate_i)
%CREATES al the times for the processing, including those to measure Peak
%Dealay and Coda Attenuation on the waveform.

f0                              =   Murat.input.originTime;
PorS                            =   Murat.input.POrS;
tCm                             =   Murat.input.startLapseTime;
vP                              =   Murat.input.averageVelocityP;
vS                              =   Murat.input.averageVelocityS;
maxtpde                         =   Murat.input.maximumPeakDelay;
tWm                             =   Murat.input.codaWindow;
    
%Conditions in case the zero time is missing in the header
if isequal(f0,[]) || isequal(eval(f0),-12345)
    
    if PorS == 2 %P-wave peak delay
        v                       =   vP;
    elseif PorS == 3 %S-wave peak delay
        v                       =   vS;
    end
    
    Distance                        =...
    sqrt((sst_i(4)-sst_i(1))^2+(sst_i(5)-sst_i(2))^2+...
    (sst_i(6)-sst_i(3))^2);

    theoreticalTime_i           =   Distance/v/1000;
    originTime_i                =   pktime_i-theoreticalTime_i;
else
    originTime_i                =   eval(f0);
end

t00                             =   tempis(1); %starting time waveform

%starting-sample of the direct-wave window
cursorPick_i                    =   floor((pktime_i-t00)*srate_i);
cursorPeakDelayMax              =   floor(cursorPick_i+maxtpde*srate_i);

%Measures for peak delay - compute shorter window if waveform is cut
if cursorPeakDelayMax > length(tempis)
    
    %set to length of seismogram
    cursorPeakDelay_i           =   length(tempis); 
    msg                         =...
        'Length of seismogram shorter than peak-delay window!';
    warning(msg);
    
else
    
    cursorPeakDelay_i           =   cursorPeakDelayMax;
    
end          

%Lapse times - in case it is not defined you are
%measuring envelopes from the peak of the direct wave, assuming that
%the entire waveform is diffusive (e.g., Deception Island).
if isequal(tCm,[])
    
    tC_i       =   (pktime_i-originTime_i)+tspm/srate_i;
    
else
    
    tC_i       =   tCm;
    
end

cursorCodaStart_i                =...
    floor((originTime_i-t00+tC_i)*srate_i-1);
cursorCodaEnd_i                  =...
    floor(cursorCodaStart_i + tWm*srate_i-1);

