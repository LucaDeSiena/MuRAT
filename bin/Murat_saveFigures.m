function  Murat_saveFigures(figureName,Path)
% function  Murat_plotSections(figure,Path)
% PLOTS sections of the desired model
%
%	Input Parameters:
%       figureName:             Matlab 3D plot
%       Path:                   name of sections (png) and figure (fig)
%
%   Output:
%       Three sections and 1 3D figure
%
SetFDefaults
saveas(figureName,Path);    
view(90,0)
ytickangle(45)
saveas(figureName,[Path '_SN'], 'png');
view(0,0)
xtickangle(45)
saveas(figureName,[Path '_WE'], 'png');
view(0,90)
saveas(figureName,[Path '_H'], 'png');
close(figureName)
end