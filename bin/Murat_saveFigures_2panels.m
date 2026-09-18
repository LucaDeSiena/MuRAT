function  Murat_saveFigures_2panels(figureName,Path)
% function  Murat_saveFigures_2panels(figureName,Path)
% PLOTS sections of the desired model in two panels
%
%	Input Parameters:
%       figureName:             Matlab 3D plot
%       Path:                   name of sections (png) and figure (fig)
%
%   Output:
%       3D figure with 2 panels
%
SetFDefaults
fig = figureName;
savefig(fig, Path);
ax1             =   subplot(1,2,1);
view(ax1,90,0)
ytickangle(45)
ax2             =   subplot(1,2,2);
view(ax2,90,0)
ytickangle(45)
saveFigureAsImage([Path '_SN']);


ax1             =   subplot(1,2,1);
view(ax1,0,0)
xtickangle(45)
ax2             =   subplot(1,2,2);
view(ax2,0,0)
xtickangle(45)
saveFigureAsImage([Path '_WE']);

ax1             =   subplot(1,2,1);
view(ax1,0,90)
ax2             =   subplot(1,2,2);
view(ax2,0,90)
saveFigureAsImage([Path '_H']);

close(figureName)

end