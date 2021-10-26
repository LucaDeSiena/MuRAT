function [image, SAChdr]    =   Murat_modv(modv,origin)
% TEST seismogram envelopes for changes in  broadening
% CREATES a figure with seismograms and envelopes for different frequencies
%
%	Input Parameters:
%       nameWaveform:           name of the SAC file
%       centralFrequencies:     vector of frequencies (Hz),if [] no filter
%       smoothingCoefficient:   coefficient to smooth envelopes
%       figOutput:              decide if you want to show figures (set 1)
%       verboseOutput:          decide if you want to show messages (set 1)
%
%   Output:
%       image:                  image with envelope at specified frequency
%       SAChdr:                 header of the SAC file
%
%   35.7 -3
% Imports SAC files
modvDEG                     =...
    sortrows([modv(:,1)-origin(1) modv(:,2)-origin(2) modv(:,3:4)],1);


modvM                       =...
    [deg2km(modvDEG(:,1))*1000,deg2km(modvDEG(:,2))*1000,modv(:,3:4)];


image                       =   [];
%% Figure
if figOutput == 1
    
    srate_i                 =   1/SAChdr.times.delta;
    
    sisma                   =   detrend(sisma,1);
    
    lsis                    =   length(sisma);
    
    tu                      =   tukeywin(lsis,0.05);
    tsisma                  =   tu.*sisma;
        
    image                   =   figure('Name',['Test Seismograms: '...
        nameWaveform],'NumberTitle','off','Position',[20,400,1200,1000]);
    lengthFrequencies       =   length(centralFrequencies);
    plotFrequencies         =   1:2:2*lengthFrequencies;
    
    if isequal(centralFrequencies,[])
        plot(times,sisma,'k-','LineWidth',2);
        xlim([SAChdr.times.a - 5 SAChdr.times.a + 20])
        SetFDefaults
    else
        for i = 1:lengthFrequencies
            % Filter creation - in loop for each frequency
            cf              =   centralFrequencies(i);
            Wn              =   ([cf-cf/3 cf+cf/3]/srate_i*2);
            [z,p,k]         =   butter(4,Wn,'bandpass');
            [sos,g]         =   zp2sos(z,p,k);
            
            fsisma          =   filtfilt(sos,g,tsisma);
            
            hsp_i           =   hilbert(fsisma);
            sp_i            =   smooth(abs(hsp_i),smoothingC/cf*srate_i);
            
            subplot(lengthFrequencies,2,plotFrequencies(i));
            plot(times,fsisma,'k-','LineWidth',2);
            xlabel('Time (s)')
            ylabel('Amplitude')
            SetFDefaults
            
            subplot(lengthFrequencies,2,plotFrequencies(i)+1);
            plot(times,sp_i,'k-','LineWidth',2);
            xlabel('Time (s)')
            ylabel('Energy')
            SetFDefaults
        end
    end
end
%%
% Checking all metadata
if verboseOutput == 1
    fprintf('<strong> Checking temporal metadata.</strong>\n');
    if isequal(SAChdr.times.o,-12345)
        disp('Origin (o) is not set.')
    else
        disp(['Origin (o) at ' num2str(SAChdr.times.o) ' s.'])
    end
    
    if isequal(SAChdr.times.a,-12345)
        disp('P wave picking (a) is not set.')
    else
        disp(['P wave picking (a) at ' num2str(SAChdr.times.a) ' s.'])
    end
    
    if isequal(SAChdr.times.t0,-12345)
        disp('S wave picking (t0) is not set.')
    else
        disp(['S wave picking (t0) at ' num2str(SAChdr.times.t0) ' s.'])
    end
    
    fprintf('<strong> Checking event location metadata.</strong>\n');
    
    if isequal(SAChdr.event.evla,-12345)
        disp('Event latitude (evla) is not set.')
    else
        disp(['Event latitude (evla) at ' num2str(SAChdr.event.evla)...
            ' degrees.'])
    end
    
    if isequal(SAChdr.event.evlo,-12345)
        disp('Event longitude (evlo) is not set.')
    else
        disp(['Event longitude (evlo) at ' num2str(SAChdr.event.evlo)...
            ' degrees.'])
    end
    
    if isequal(SAChdr.event.evdp,-12345)
        disp('Event depth (evdp) is not set.')
    else
        disp(['Event depth (evdp) at ' num2str(SAChdr.event.evdp) ' km.'])
    end
    
    
    fprintf('<strong> Checking station location metadata.</strong>\n');
    
    if isequal(SAChdr.station.stla,-12345)
        disp('Station latitude (stla) is not set.')
    else
        disp(['Station latitude (stla) at ' num2str(SAChdr.station.stla)...
            ' degrees.'])
    end
    
    if isequal(SAChdr.station.stlo,-12345)
        disp('Station longitude (stlo) is not set.')
    else
        disp(['Station longitude (stlo) at ' num2str(SAChdr.station.stlo)...
            ' degrees.'])
    end
    
    if isequal(SAChdr.station.stel,-12345)
        disp('Station elevation (stel) is not set.')
    else
        disp(['Station elevation (stel) at ' num2str(SAChdr.station.stel)...
            ' m.'])
    end
end
function SetFDefaults()
% DEFAULT settings for MuRAt figures
ax = gca;
ax.GridColor = [0 0 0];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.3;
ax.LineWidth = 1.5;
ax.FontSize  = 12;
grid on