function [sp_i,hsp_i,SAChdr_i,srate_i]...
    =   Murat_envelope(Murat,listSac_i)
%ENVELOPE calculation, the code also outputs the picking time 

cf                      =   Murat.input.centralFrequency;
nf                      =   Murat.input.smoothing;

% Imports SAC files
[~,sisma,SAChdr]        =   fget_sac(listSac_i); 
% Sampling rate
srate_i                 =   1/SAChdr.times.delta; %sampling frequency
% Remove mean and trends
sisma                   =   detrend(sisma,1);%remove trend
% Seismogram length
lsis                    =   length(sisma);
% Tapering
tu                      =   tukeywin(lsis,0.05);%avoid windowing
tsisma                  =   tu.*sisma;

%% Calculations

% Filter creation - in loop as sampling might change
Wn                      =   ([cf-cf/3 cf+cf/3]/srate_i*2); %frequency band
[z,p,k]                 =   butter(4,Wn,'bandpass'); %butter filter
[sos,g]                 =   zp2sos(z,p,k); % Convert to SOS form

fsisma                  =   filtfilt(sos,g,tsisma);% filtered waveform

SAChdr_i                =   SAChdr;
hsp_i                   =   hilbert(fsisma); %saved for CN method
sp_i                    =   smooth(abs(hsp_i),nf/cf*srate_i);%envelope
