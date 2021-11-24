% function  Murat_plotAfter(Murat)
% PLOTS sections of the desired model with specified 
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
[file,path]                 =   uigetfile('*.fig');

if isequal(file,0)
   error('User selected Cancel!');   
else
   disp(['User selected ', fullfile(path,file)]);
end
figure                      =   openfig(fullfile(path,file),'visible');
SetFDefaults
limits                      =   input('Axes limits ([min max])?');
caxis(limits)

function SetFDefaults()
% DEFAULT settings for MuRAt figures
ax = gca;
ax.GridColor = [0 0 0];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.3;
ax.LineWidth = 1.5;
ax.FontSize  = 12;
grid on
end