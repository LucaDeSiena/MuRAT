function [cursorPick_i, pktime_i, v_i]  =...
    Murat_picking(tempis,PTime,STime,PorS,vP,vS,srate_i,SAChdr) %#ok<INUSD>
% CHECKS if you are working with P or S picking and if this picking is
% inside the waveform.
%
% Input parameters:
%    tempis:        times from seismograms
%    PTime:         P-wave time on the waveform
%    STime:         S-wave time on the waveform
%    PorS:          As defined by the used
%    vP:            P-wave velocity 
%    vS:            S-wave velocity model
%    srate_i:       sampling rate
%    SAChdr:        SAC header from trace rate
%
% Output parameters:
%    cursorPick_i:  position of the picking on the trace
%    pktime_i:      picking in seconds
%    v_i:           chosen average velocity

if PorS == 2
    pktime_i                            =   eval(PTime);
    v_i                                 =   vP;
elseif PorS == 3
    pktime_i                            =   eval(STime);
    v_i                                 =   vS;
end

% If picking is not on the waveform
if pktime_i < tempis(1)
    error(['The picking is set before the start of the recording number '...
        num2str(i_label)])
elseif pktime_i > tempis(end)
    error(['The picking is set after the end of the recording number '...
        num2str(i_label)])
end

t00                                     =   tempis(1);

cursorPick_i                            =   floor((pktime_i-t00)*srate_i);
end



