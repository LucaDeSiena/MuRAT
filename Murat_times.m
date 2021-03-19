function [theoreticalTime_i, tC_i, cursorPick_i, cursorPeakDelay_i,...
    cursorCodaStart_i, cursorCodaEnd_i] =...
    Murat_times(Murat,sst_i,listSac_i)
%CREATES al the times for the processing, including those to measure Peak
%Dealay and Coda Attenuation on the waveform.

f0                              =   Murat.input.originTime;
fP                              =   Murat.input.PTime;
fS                              =   Murat.input.STime;
PorS                            =   Murat.input.POrS;
tCm                             =   Murat.input.startLapseTime;
vP                              =   Murat.input.averageVelocityP;
vS                              =   Murat.input.averageVelocityS;
maxtpde                         =   Murat.input.maximumPeakDelay;
tWm                             =   Murat.input.codaWindow;
 
% Imports SAC files
[tempis,~,SAChdr]               =   fget_sac(listSac_i); 

% Sampling rate
srate_i                         =   1/SAChdr.times.delta;

%Check if you are doing P or S peak delay
if PorS == 2 %P-wave peak delay
    pktime_i                    =   eval(fP);%improve here
    v                           =   vP;
elseif PorS == 3 %S-wave peak delay
    pktime_i                    =   eval(fS);%improve here
    v                           =   vS;
end

%If picking is not on the waveform
if pktime_i < tempis(1)
    error('The picking is set before the start of the recording!')
elseif pktime_i > tempis(end)
    error('The picking is set after the end of the recording!')
end


%Conditions in case the zero time is missing in the header
if isequal(f0,[]) || isequal(eval(f0),-12345)
    
    Distance                    =...
    sqrt((sst_i(4)-sst_i(1))^2+(sst_i(5)-sst_i(2))^2+...
    (sst_i(6)-sst_i(3))^2);

    theoreticalTime_i           =   Distance/v/1000;
    originTime_i                =   pktime_i-theoreticalTime_i;
    
else
    
    originTime_i                =   eval(f0);
    theoreticalTime_i           =   pktime_i - originTime_i;
    
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
    
    tC_i                        =   (pktime_i-originTime_i)+tspm/srate_i;
    
else
    
    tC_i                        =   tCm;
    
end

cursorCodaStart_i               =...
    floor((originTime_i-t00+tC_i)*srate_i-1);
cursorCodaEnd_i                 =...
    floor(cursorCodaStart_i + tWm*srate_i-1);

%Setting the end of the coda window as the end of the
%waveform (in case the seismogram is too short)

if cursorCodaEnd_i > length(tempis)
    cursorCodaEnd_i             = length(tempis);
end
