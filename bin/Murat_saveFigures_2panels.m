function  Murat_saveFigures_2panels(figureName,Path)
% function  Murat_plotSections(figure,Path)
% PLOTS sections of the desired model
%
%	Input Parameters:
%       figureName:             Matlab 3D plot
%       Path:                   name of sections (png) and figure (fig)
%
%   Output:
%       Three sections and 1 3D figure with 2 panels
%
SetFDefaults
saveas(figureName,Path);    

ax1             =   subplot(1,2,1);
view(ax1,90,0)
ytickangle(45)
ax2             =   subplot(1,2,2);
view(ax2,90,0)
ytickangle(45)
saveas(figureName,[Path '_SN'], 'png');

ax1             =   subplot(1,2,1);
view(ax1,0,0)
xtickangle(45)
ax2             =   subplot(1,2,2);
view(ax2,0,0)
xtickangle(45)
saveas(figureName,[Path '_WE'], 'png');

ax1             =   subplot(1,2,1);
view(ax1,0,90)
ax2             =   subplot(1,2,2);
view(ax2,0,90)
saveas(figureName,[Path '_H'], 'png');

close(figureName)

end