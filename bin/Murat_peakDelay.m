function peakDelay_i        =...
    Murat_peakDelay(sp_i,cursorPick_i,srate_i,cursorPeakDelay_i)
% function peakDelay_i            =...
%     Murat_peakDelay(sp_i,cursorPick_i,srate_i,cursorPeakDelay_i)
%
% CALCULATES peak delay time
%
% Input parameters:
%    sp_i:              envelope at all frequencies
%    cursorPick_i:      pick time for the chosen phase on trace
%    srate_i:           sampling rate
%    cursorPeakDelay_i: maximum pick delay on trace
%
% Output parameters:
%    peakDelay_i:       peak delay on trace

lsp                         =   size(sp_i,2);
peakDelay_i                 =   zeros(lsp,1);

for i = 1:lsp
    pdsp                    =   sp_i(cursorPick_i:cursorPeakDelay_i,i);

    [~,timeMaxAmplitude]    =	max(pdsp);

    peakDelay_i(i)          =   timeMaxAmplitude/srate_i;
end

end

