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
[file,path]                 =   uigetfile('*.mat');

if isequal(file,0)
   error('User selected Cancel!');   
else
   disp(['User selected ', fullfile(path,file)]);
end
load(fullfile(path,file));
