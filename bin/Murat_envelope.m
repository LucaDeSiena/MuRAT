function [tempis, sp_i] = Murat_envelope(cf, listSac_i, envelopeSmoothTime)
% MURAT_ENVELOPE
% Compute the signal envelope for all frequency bands.
%
% Input parameters:
%    cf:                  central frequencies
%    listSac_i:           SAC file name
%    envelopeSmoothTime:  envelope smoothing window in seconds
%
% Output parameters:
%    tempis:              seismogram time vector
%    sp_i:                squared envelopes for the different frequencies

%%% MODIFICATION:
%%% Added envelopeSmoothTime as a third input argument.
%%% If not provided, a default value of 1.0 s is used
%%% to preserve the original behavior.

if nargin < 3 || isempty(envelopeSmoothTime)
    envelopeSmoothTime = 1.0;
end

lcf = length(cf);

[tempis, sisma, SAChdr_i] = fget_sac(listSac_i);
srate_i = 1 / SAChdr_i.times.delta;

sisma = detrend(sisma,1);
lsis  = length(sisma);

tu     = tukeywin(lsis, 0.05);
tsisma = tu .* sisma;

sp_i = zeros(lsis, lcf);

%%% MODIFICATION:
%%% Improved validation of the sampling rate.
%%% Previously, only the value -12345 was checked.
if isequal(srate_i, -12345) || ~isfinite(srate_i) || srate_i <= 0
    error(['Waveform ' listSac_i ' has no valid sampling rate!'])
end

%%% MODIFICATION:
%%% Convert the smoothing time (in seconds) to a number of samples.
%%% Example: 0.25 s * 100 Hz = 25 samples.
winSamples = max(1, round(envelopeSmoothTime * srate_i));

for i = 1:lcf

    % Create the bandpass filter for the current frequency
    Wn = [cf(i)-cf(i)/3, cf(i)+cf(i)/3] / srate_i * 2;

    %%% MODIFICATION:
    %%% Ensure a valid frequency band to avoid butter filter errors.
    if Wn(1) <= 0
        Wn(1) = 0.001;
    end
    if Wn(2) >= 1
        Wn(2) = 0.999;
    end
    if Wn(2) <= Wn(1)
        warning('Skipping frequency %.2f Hz for %s: invalid bandpass.', cf(i), listSac_i);
        continue
    end

    [z,p,k] = butter(4, Wn, 'bandpass');
    [sos,g] = zp2sos(z,p,k);

    fsisma = filtfilt(sos, g, tsisma);

    %%% MODIFICATION:
    %%% Previously, MuRAT used round(srate_i), corresponding
    %%% to an approximately 1-second smoothing window.
    %%% The window length is now controlled by the input parameter.
    [sp, ~] = envelope(fsisma, winSamples, 'rms');

    % Compute energy as squared envelope
    sp_i(:,i) = sp.^2;
end

end