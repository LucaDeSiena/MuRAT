function Murat_hitmap(CoordGrid,rayCrossing,name)
% 
% HITMAP produces horizontal and vertical cross section displaying the
%                   seismic ray count as hitmap. 3 default sections in the
%                   model are built (N-S, W-E, and horizontal). 
%
% Input parameters:
%    CoordGrid:     3D coordinates in matrix format
%    rayCrossing:   vector containing number of rays per cell
%
x = CoordGrid(:,1);
y = CoordGrid(:,2);
z = CoordGrid(:,3);
tbl = table(x,y,z,rayCrossing);

% Build default hitmaps
myfig();
h = heatmap(tbl,'x','y','ColorVariable','rayCrossing','Colormap', hot);
h.NodeChildren(3).YDir='normal';
saveFigureAsImage([name '_H']);
close
myfig();
vx = heatmap(tbl,'x','z','ColorVariable','rayCrossing','Colormap', hot);
vx.NodeChildren(3).YDir='normal';
saveFigureAsImage([name '_WE']);
close
myfig();
vy = heatmap(tbl,'y','z','ColorVariable','rayCrossing','Colormap', hot);
vy.NodeChildren(3).YDir='normal';
saveFigureAsImage([name '_SN']);
close