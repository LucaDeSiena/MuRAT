function peakDelay_i    =...
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

% Extract the windowed matrix: rows from pick to peakDelay, all columns
window      =   sp_i(cursorPick_i:cursorPeakDelay_i, :);

% Find index of max within each column (relative to window start)
[~, idx]    =   max(window, [], 1);

% Convert to seconds and return as column vector
peakDelay_i =   (idx.' ) / srate_i;

end

