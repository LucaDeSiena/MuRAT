function  Murat             =   Murat_plotMore(sections)
% function  Murat             =   Murat_plotAfter(Murat)
% PLOTS sections of the desired model with specified sections
%
%
[file,path]                 =   uigetfile('*.mat');

if isequal(file,0)
   error('User selected Cancel!');   
else
   disp(['User selected ', fullfile(path,file)]);
end
Murat                       =   load(fullfile(path,file));
Murat.input.sections        =   sections;
Murat                       =   Murat_plot(Murat);
