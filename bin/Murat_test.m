function [image, SAChdr]    =...
    Murat_test(nameWaveform,centralFrequencies,smoothingCoefficient)
% TEST seismogram envelopes for changes in  broadening
% CREATES a figure with seismograms and envelopes for different frequencies
%
%	Input Parameters:
%       nameWaveform:           name of the SAC file
%       centralFrequencies:     vector of frequencies (Hz),if [] no filter
%       smoothingCoefficient:   coefficient to smooth envelopes
%
%   Output:
%       image:                  image with envelope at specified frequency
%

% Imports SAC files
[times,sisma,SAChdr]        =   fget_sac(nameWaveform);
% Sampling rate
srate_i                     =   1/SAChdr.times.delta;
% Remove mean and trends
sisma                       =   detrend(sisma,1);
% Seismogram length
lsis                        =   length(sisma);
% Tapering
tu                          =   tukeywin(lsis,0.05);
tsisma                      =   tu.*sisma;

%% Figure
image                       =   figure('Name',['Test Seismograms: '...
    nameWaveform],'NumberTitle','off','Position',[20,400,1200,1000]);
lengthFrequencies           =   length(centralFrequencies);
plotFrequencies             =   1:2:2*lengthFrequencies;

if isequal(centralFrequencies,[])
    plot(times,sisma);
    xlim([SAChdr.times.a - 5 SAChdr.times.a + 20])
    SetFDefaults
else
    for i = 1:lengthFrequencies
        % Filter creation - in loop for each frequency
        cf                      =   centralFrequencies(i);
        Wn                      =   ([cf-cf/3 cf+cf/3]/srate_i*2); %frequency band
        [z,p,k]                 =   butter(4,Wn,'bandpass'); %butter filter
        [sos,g]                 =   zp2sos(z,p,k); % Convert to SOS form
        
        fsisma                  =   filtfilt(sos,g,tsisma);% filtered waveform
        
        hsp_i                   =   hilbert(fsisma); %saved for CN method
        sp_i                    =...
            smooth(abs(hsp_i),smoothingCoefficient/cf*srate_i);%envelope
        
        subplot(lengthFrequencies,2,plotFrequencies(i));
        plot(times,fsisma);
        xlabel('Time (s)')
        ylabel('Amplitude')
        SetFDefaults
        
        subplot(lengthFrequencies,2,plotFrequencies(i)+1);
        plot(times,sp_i);
        xlabel('Time (s)')
        ylabel('Energy')
        SetFDefaults
    end
end

if isequal(SAChdr.times.o,-12345)
    disp('Origin (o) is not set.')
else
    disp(['Origin (o) at ' num2str(SAChdr.times.o) ' s.'])
end

if isequal(SAChdr.times.a,-12345)
    disp(['P wave picking (a) is not set.'])
else
    disp(['P wave picking (a) at ' num2str(SAChdr.times.a) ' s.'])
end

if isequal(SAChdr.times.t0,-12345)
    disp('S wave picking (t0) is not set.')
else
    disp(['S wave picking (t0) at ' num2str(SAChdr.times.t0) ' s.'])
end


