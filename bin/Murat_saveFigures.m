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
SetFDefaults()
savefig(figureName,Path);    
view(90,0)
ytickangle(45)
saveFigureAsImage([Path '_SN']);
view(0,0)
xtickangle(45)
saveFigureAsImage([Path '_WE']);
view(0,90)
saveFigureAsImage([Path '_H']);

close(figureName)
end