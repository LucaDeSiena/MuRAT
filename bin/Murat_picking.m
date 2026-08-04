function [cursorPick_i, pktime_i, v_i]  =   Murat_picking(tempis,PTime,...
    STime,PorS,vP,vS,srate_i,listaSAC_i,SAChdr) %#ok<INUSD>
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
%
% Output parameters:
%    cursorPick_i:  position of the picking on the trace
%    pktime_i:      picking in seconds
%    v_i:           chosen average velocity

if PorS == 2
    pktime_i    =   eval(PTime);
    v_i         =   vP;
elseif PorS == 3
    pktime_i    =   eval(STime);
    v_i         =   vS;
end

% If picking is not on the waveform
if pktime_i < tempis(1)
    error(['Picking is before the start of the waveform ' listaSAC_i])
elseif pktime_i > tempis(end)
    error(['Picking is after the end of the waveform ' listaSAC_i])
end

cursorPick_i    =   floor((pktime_i-tempis(1))*srate_i);
end